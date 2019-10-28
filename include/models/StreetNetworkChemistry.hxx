#ifndef POLYPHEMUS_FILE_MODELS_STREETNETWORKCHEMISTRY_HXX

#include "StreetNetworkTransport.cxx"
#include "StreetTransport.cxx"

namespace Polyphemus
{


  using namespace std;
  using namespace blitz;


  //////////////////
  // STREET //
  //////////////////


  //! StreetNetworkChemistry model.
  template<class T, class ClassChemistry>
  class StreetNetworkChemistry: public StreetNetworkTransport<T>
  {

  protected:

    static const T pi;

    /*** Meteorological data ***/
    Array<T, 2> attenuation;
    Array<T, 2> specific_humidity;
    Array<T, 2> pressure;
    Array<T, 2> temperature;

    //***** Chemistry *****//

    ClassChemistry Chemistry_;
    string option_chemistry;

    //*** Photolysis ***//
    //! Number of photolysis reactions.
    int Nr_photolysis;
    //! Grid for photolysis reactions.
    RegularGrid<T> GridR_photolysis;
    //! List of species with photolysis reactions.
    vector<string> photolysis_reaction_list;
    //! List of photolysis rates.
    Array<T, 1> photolysis_rate_;
    //! List of altitudes at which photolysis rates are available.
    vector<string> altitudes_photolysis;

    //! Starting date of photolysis-rates files.
    Date photolysis_date_min;
    //! Time step of photolysis-rates files in seconds.
    T photolysis_delta_t;
    //! Number of days in photolysis-rates files.
    int Nphotolysis_days;

    //! First time angle in photolysis-rates data.
    T photolysis_time_angle_min;
    //! Time-angle step in photolysis-rates data.
    T photolysis_delta_time_angle;
    //! Number of time angles in photolysis-rates data.
    int Nphotolysis_time_angle;
    //! Grid for time angles of photolysis rates.
    RegularGrid<T> Grid_time_angle_photolysis;

    //! First latitude in photolysis-rates data.
    T photolysis_latitude_min;
    //! Latitude step in photolysis-rates data.
    T photolysis_delta_latitude;
    //! Number of latitudes in photolysis-rates data.
    int Nphotolysis_latitude;
    //! Grid for latitudes of photolysis rates.
    RegularGrid<T> Grid_latitude_photolysis;

    //! Number of vertical levels in photolysis-rates data.
    int Nphotolysis_z;
    //! Grid for altitudes of photolysis rates.
    RegularGrid<T> GridZ_photolysis;

  public:

    /*** Constructors and destructor ***/

    StreetNetworkChemistry();
    StreetNetworkChemistry(string config_file);
    virtual ~StreetNetworkChemistry();

    /*** Initializations ***/

    void ReadConfiguration();
    void CheckConfiguration();
    void Allocate();    
    void InitPhotolysis(Date date);

    void Init();
    void InitStreet();
    void InitStep();
    void InitData();
    void InitData(string input_file, Array<T, 2>& input_data);

    /*** Access Methods ***/

    void SetStreetAdditionalMeteo(int street_index, T attenuation,
                                  T specific_humidity, T pressure,
                                  T temperature);

    /*** Computational Methods ***/
    
    void Forward();

    //LL: Remove stationary regime*****************************************
    void ComputeStreetConcentrationNoStationary();
    void Chemistry(Date current_date_tmp,
		   Array<T, 1>& concentration_array,
		   const T attenuation_,
		   const T specific_humidity_,
		   const T pressure_,
		   const T temperature_,
		   const T longitude_,
		   const T latitude_,
		   const Array<T, 1> photolysis_rate,
		   const T sub_delta_t);
    //***********************************************************************


    /*** Access Methods: Chemistry ***/
    void Chemistry();
    void InitChemistry();
    int GetNr_photolysis() const {return Nr_photolysis;}
    vector<string> GetPhotolysisReactionList() const 
    {return photolysis_reaction_list;}
    int GetNs_source() const {return 0;}
    int GetNz_source() const {return 0;}
    int SourceGlobalIndex(int s) const {return 0;}

    bool WithChemistry();
    string GetChemicalMechanism();
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKCHEMISTRY_HXX
#endif
