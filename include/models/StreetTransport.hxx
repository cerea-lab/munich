#ifndef POLYPHEMUS_FILE_MODELS_STREETTRANSPORT_HXX

#include <list>

namespace Polyphemus
{

  //////////////
  // INCLUDES //
  //////////////
  
  using namespace std;

  ////////////
  // STREET //
  ////////////


  //! Class that stores all data about a street.
  template<class T>
  class Street
  {

  protected:

    //! Street id
    int street_id_;

    //! Longitude of the street center.
    T longitude_;

    //! Latitude of the street center.
    T latitude_;

    //! Begin intersection id
    int begin_inter_;

    //! End intersection id
    int end_inter_;
    
    //! Street length (m)
    T length_;

    //! Street width (m)
    T width_;

    //! Building height in the street (m)
    T height_;

    //! Street angle (rad). 
    T street_angle_;

    //! Street angle for the intersection (rad)
    T street_angle_inter_;

    //! Street wind velocity (m/s)
    T ustreet_;

    //! Outgoing volume flux from the street (m3/s).
    T outgoing_flux_;

    //! Incoming volume flux from other streets (m3/s).
    T incoming_flux_;

    //! Deposition rate (ug/s)
    T deposition_rate_;

    //! Turbulent transfer velocity (m/s)
    T transfer_velocity_;

    //! Standard deviation of vertical wind speed (m/s)
    T sigma_w_;

    //! Number of species.
    int ns_local_;

    //! Emission rate (ug/s). 
    Array<T, 1> emission_;

    //! Inflow rate to the street (ug/s).
    Array<T, 1> inflow_rate_;

    //! Outflow rate to the atmosphere (ug/s)
    Array<T, 1> massflux_roof_to_background_;

    //! Inflow rate from the atmosphere to the canyon through the roof (ug/s)
    Array<T, 1> massflux_background_to_roof_;

    //! Outflow rate to the atmosphere through the intersection (ug/s)
    Array<T, 1> massflux_to_background_;

    //! Inflow rate from the atmosphere to the canyon through the intersection (ug/s).
    Array<T, 1> massflux_from_background_;

    //! List of species concentrations (ug/m3).
    Array<T, 1> concentration_;

    //! List of species initial concentrations (ug/m3).
    Array<T, 1> initial_concentration_;

    //! List of species mass change (ug).
    Array<T, 1> street_quantity_delta_;

    //! List of background species concentrations (ug/m3).
    Array<T, 1> background_concentration_;
   
    //! If all concentrations in the street reach the stationary state.
    bool is_stationary_;

    //! If the local wind speed is lower than the threshold.
    bool is_low_wind_;

    //! Meteo data in the background grid cell.
    T wind_direction_;
    T wind_speed_;
    T pblh_;
    T ust_;
    T lmo_;

    //! Meteo data for chemisty.
    T attenuation_;
    T specific_humidity_;
    T pressure_;
    T temperature_;

    //! Photolysis rate (1/s)
    Array<T, 1> photolysis_rate_;

  public:

     /*** Constructor and destructor ***/

    Street(int street_id, int begin_inter, int end_inter, 
           T length, T width, T height, int ns_local, int nr_photolysis);
    virtual ~Street();   

    /*** Methods ***/

    int GetStreetID() const;
    int GetBeginIntersectionID() const;
    int GetEndIntersectionID() const;
    T GetLongitude() const;
    T GetLatitude() const;
    void SetCoordinate(T longitude, T latitude);
    T GetPhotolysisRate(int r) const;
    void SetPhotolysisRate(Array<T, 1> photolysis_rate);
    T GetLength() const;
    T GetWidth() const;
    T GetHeight() const;
    T GetStreetAngle() const;
    void SetStreetAngle(T street_angle);
    T GetStreetAngleIntersection() const;
    void SetStreetAngleIntersection(T street_angle);
    T GetStreetWindSpeed() const;
    void SetStreetWindSpeed(T ustreet);
    T GetStreetUstar() const;
    void SetStreetUstar(T ustar_street);
    T GetStreetVolume() const;
    T GetStreetQuantity(int s) const;
    T GetStreetConcentration(int s) const;
    void SetStreetConcentration(T concentration, int s);
    T GetInitialStreetConcentration(int s) const;
    void SetInitialStreetConcentration(T concentration, int s);
    T GetStreetQuantityDelta(int s) const;
    void SetStreetQuantityDelta(T street_quantity_delta, int s);
    T GetBackgroundConcentration(int s) const;
    void SetBackgroundConcentration(T background_concentration, int s);
    T GetEmission(int s) const;
    void SetEmission(Array<T, 1> emission);
    T GetWindDirection() const;
    T GetPBLH() const;
    T GetLMO() const;
    void SetMeteo(T wind_direction, T wind_speed, 
                  T pblh, T ust, T lmo);
    T GetAttenuation() const;
    T GetSpecificHumidity() const;
    T GetPressure() const;
    T GetTemperature() const;
    void SetMeteoChemistry(T attenuation, T specific_humidity, 
                  T pressure, T temperature);
    T GetInflowRate(int s) const;
    void SetInflowRate(T inflow_rate, int s);
    T GetOutgoingFlux() const;
    void SetOutgoingFlux(T outgoing_flux);
    T GetIncomingFlux() const;
    void SetIncomingFlux(T incoming_flux);
    T GetMassfluxRoofToBackground(int s) const;
    void SetMassfluxRoofToBackground(T massflux_roof_to_background, int s);
    T GetMassfluxBackgroundToRoof(int s) const;
    void SetMassfluxBackgroundToRoof(T massflux_background_to_roof, int s);
    T GetMassfluxToBackground(int s) const;
    void SetMassfluxToBackground(T massflux_to_background, int s);
    T GetMassfluxFromBackground(int s) const;
    void SetMassfluxFromBackground(T massflux_from_background, int s);
    T GetDepositionRate() const;
    void SetDepositionRate(T deposition_rate);
    T GetTransferVelocity() const;
    void SetTransferVelocity(T transfer_velocity);
    T GetSigmaW() const;
    void SetSigmaW(T sigma_w);
    bool GetStationary() const;
    void SetStationary(bool is_stationary);

  };


  //! Class that stores all data about an intersection.
  template<class T>
  class Intersection
  {

  protected:

    //! Intersection id
    int id_;

    //! Longitude
    T x_;

    //! Latitude
    T y_;

    //! Number of streets which are connected to the intersection
    int nstreet_;

    //! List of the streets
    Array<int, 1> street_list_;

    //! If the intersection is artificially created. It is created because the street cross two grid cells.
    bool is_virtual_;

    //! Flux matrix
    Array<T, 2> flux_matrix_;

    //! Flux matrix using the Gaussian distribution
    Array<T, 2> gaussian_matrix_;

    //! Meteo data in the background grid cell.
    T wind_direction_;
    T wind_speed_;
    T pblh_;
    T ust_;
    T lmo_;
    T sigma_theta_;
    Array<T, 1> theta_, gaussian_;

  public:

     /*** Constructor and destructor ***/

    Intersection(int id, T x, T y, int nstreet,
                 Array<int, 1> street_list, bool is_virtual,
                 Array<T, 2> flux_matrix, Array<T, 2> gaussian_matrix);
    virtual ~Intersection();   

    /*** Methods ***/

    int GetID() const;
    T GetX() const;
    T GetY() const;
    T GetNStreet() const;
    Array<int, 1>& GetStreetList();
    bool IsVirtual();
    Array<T, 2>& GetFluxMatrix();
    void SetFluxMatrix(Array<T, 2> flux_matrix);
    Array<T, 2>& GetGaussianMatrix();
    void SetGaussianMatrix(Array<T, 2> gaussian_matrix);
    void SetMeteo(T wind_direction, T wind_speed, 
                  T pblh, T ust, T lmo);    
    T GetPBLH() const;
    T GetUST() const;
    T GetLMO() const;
    T GetWindSpeed() const;
    T GetWindDirection() const;
    T GetSigmaTheta() const;
    void SetSigmaTheta(T sigma_theta);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETTRANSPORT_HXX
#endif
