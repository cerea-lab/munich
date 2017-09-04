#ifndef POLYPHEMUS_FILE_MODELS_STREETNETWORKTRANSPORT_HXX

#include "BaseModel.cxx"
#include "StreetTransport.cxx"

namespace Polyphemus
{


  using namespace std;
  using namespace blitz;


  ////////////////////////////
  // StreetNetworkTransport //
  ////////////////////////////


  /*! StreetNetworkTransport model describes the atmospheric concentrations of 
    pollutants in an urban street network.
   */
  template<class T>
  class StreetNetworkTransport: public BaseModel<T>
  {

  protected:

    static const T pi, earth_radius, karman;

    /*** Configurations ***/

    string option_transfer, option_ustreet;
    bool is_stationary;
    //! Output configuration.
    string output_config;
    string output_dir;
    bool text_file;

    /*** Domain ***/
    //! 2D grid for streets.
    RegularGrid<T> GridST2D;

    //! 2D grid for species.
    RegularGrid<T> GridS2D;

    /*** Emission ***/

    int Nt_emis;
    int Ns_emis;
    //! Emission rate
    Array<T, 3> emission;
    //! List of species with emissions.
    vector<string> species_list_emis;

    /*** Street data ***/

    string file_street;
    int total_nstreet;
    Array<int, 1> id_street, begin_inter, end_inter;
    Array<T, 1> length, width, height;
    //! Pointer to the current street.
    typename vector<Street<T>* >::iterator current_street;
    //! Street list.
    vector<Street<T>*> StreetVector, StreetVectorInter;

    /*** Background concentration data ***/

    int Nt_background;
    int Ns_background;
    Array<T, 3> background_concentration;
    //! List of species with background concentrations.
    vector<string> species_list_background;
    T cell_volume;

    /*** Meteorological data ***/

    int Nt_meteo; 
    T ustreet_min; // Minimum wind speed in the streets (m/s)
    Array<T, 2> wind_direction_arr; // Wind direction (rad)
    Array<T, 2> wind_speed_arr; // Wind speed  (m/s)
    Array<T, 2> pblh_arr; // Planetary boundary layer height (m)
    Array<T, 2> ust_arr; // Friction velocity (m/s)
    Array<T, 2> lmo_arr; // Monin-Obukhov length (m)
    int meteo_index;

    /*** Intersection data ***/

    string file_intersection;
    //! Number of intersections.
    int nintersection;
    //! ID of the intersection.
    Array<int, 1> id_inter;
    //! Number of streets which are connected to the intersection.
    Array<int, 1> nstreet;
    //! Maximum number of streets which can be connected to the intersection.
    int maxnstreet;
    Array<T, 1> x_inter, y_inter;
    //! List of the streets which are connected to the intersection.
    Array<int, 2> street_list; 
    Array<bool, 1> is_virtual;
    Array<T, 2> wind_direction_inter, wind_speed_inter, pblh_inter, ust_inter, lmo_inter;
    //! Pointer to the current intersection.
    typename vector<Intersection<T>* >::iterator current_intersection;
    //! Intersections list.
    vector<Intersection<T>*> IntersectionVector;

    /*** Concentration ***/
    Data<T, 2> StreetConcentration;

  public:

    /*** Constructors and destructor ***/

    StreetNetworkTransport();
    StreetNetworkTransport(string config_file);
    virtual ~StreetNetworkTransport();

    /*** Initializations ***/

    void ReadConfiguration();
    void ReadStreetData();
    void CheckConfiguration();
    
    void Init();
    void InitStep();
    void InitData(string input_file, Array<T, 2>& input_data);
    void InitData();

    void InitStreet();
    void EraseStreet();
    void ClearStreetVector();

    void InitIntersection();
    void EraseIntersection();
    void ClearIntersectionVector();

    void Allocate();

    /*** Access Methods ***/

    void SetStreetConcentration();
    Data<T, 2>& GetStreetConcentration();
    void SetInitialStreetConcentration();
    void SetStreetBackgroundConcentration(int street_index,
                                          Array<T, 1> background_concentration);
    int GetNStreet() const;
    int GetNumberIntersection() const;
    void GetStreetCoordinate(int street_index, T& longitude, T& latitude);
    void SetStreetMeteo(int street_index, T wind_direction, T wind_speed, 
                        T pblh, T ust, T lmo);
    void SetDelta_t(T delta_t);
    T GetDelta_t() const;
    void SetDateMin(Date date_min);
    Date GetDateMin() const;
    void GetIntersectionCoordinate(int intersection_index, T& longitude, T& latitude);
    void SetIntersectionMeteo(int intersection_index, T wind_direction, T wind_speed, 
                              T pblh, T ust, T lmo);
    void SetCurrentStreet(int index);
    void SetCurrentIntersection(int index);
    T GetStreetQuantity(int index, int s);
    T GetStreetHeight(int index);
    T GetStreetLength(int index);
    int GetStreetID(int index);
    T GetStreetVolume(int index);
    T GetMassTransferBackground(int index, int s);
    T GetMassFluxExchange(int index, int s);
    T GetStreetEmissionRate(int index, int s);

    int GetIntersectionID(int index);

    /*** Computational Methods ***/
    
    void Forward();
    void Transport();

    void ComputeIntersection(T wind_dir_inter, Intersection<T>* intersection);
    void ComputeIntersectionFlux(Array<T, 2>& extended_matrix,
                                 T wind_dir_inter);
    void ComputeGaussianFluxMatrix(T gaussian_factor, Intersection<T>* intersection);
    void CreateExtendedMatrix(Array<int, 1> ind_street_in,
                              Array<int, 1> ind_street_out,
                              Array<T, 2> flux_matrix,
                              Array<T, 2>& extended_matrix);
    void ComputeStreetAngle();
    T ComputeSigmaV(T lmo, T pblh, T ustar);
    T ComputeGaussian(double theta, double sigma_theta, 
                      double theta0);
    void ComputeSigmaW();
    void ComputeUstreet();
    void ComputeWindDirectionFluctuation();
    T ComputeUstreetSIRANE();
    T ComputeUstreetLemonsu();
    T ComputeUH(T c, T ustar_street);
    T ComputeBesselC(T z0_build, T delta_i);
    void GetMin(Array<T, 1> arr, int length, T& minimum, int& index);
    void ComputeTransferVelocity();
    void ComputeAlpha(int nstreet_in, int nstreet_out, 
                      Array<double, 1> P_in, Array<double, 1> P_out,
                      Array<double, 2>& alpha, Array<double, 2>& flux_matrix);
    void ComputeInflowRateExtended();
    void InitInflowRate();
    void ComputeBackgroundConcentration();
    void ComputeStreetConcentration();
    void IsStationary(bool& is_stationary);
    void OutputSaver();
    void InitOutputSaver();
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKTRANSPORT_HXX
#endif
