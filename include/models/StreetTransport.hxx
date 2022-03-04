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

    //! Building height in the street (m)
    int typo_;
    
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
    T rain_;
    T liquidwatercontent_;
    T cloud_height_;
    T pH_;
    T pressure_;
    T temperature_;
    T relative_humidity_;
    T specific_humidity_;
    
    /*** Dry deposition ***/
    Array<T, 1> street_dry_deposition_velocity_;
    Array<T, 1> street_dry_deposition_flux_;
    Array<T, 1> wall_dry_deposition_velocity_;
    Array<T, 1> wall_dry_deposition_flux_;
    Array<T, 1> roof_dry_deposition_velocity_;

    /*** Scavenging ***/
    Array<T, 1> street_scavenging_flux_overcanopy_;
    Array<T, 1> street_scavenging_coefficient_;
    Array<T, 1> street_scavenging_flux_;    

  public:

     /*** Constructor and destructor ***/

    // Street(int street_id, int begin_inter, int end_inter, 
    //        T length, T width, T height, int ns_local);
    Street(int street_id,
	   int begin_inter,
	   int end_inter, 
           T length,
	   T width,
	   T height,
	   int typo,
	   int ns_local);
    
    virtual ~Street();   

    /*** Methods ***/

    int GetStreetID() const;
    int GetBeginIntersectionID() const;
    int GetEndIntersectionID() const;
    T GetLongitude() const;
    T GetLatitude() const;
    void SetCoordinate(T longitude, T latitude);
    virtual T GetPhotolysisRate(int r) const;
    virtual void SetPhotolysisRate(Array<T, 1> photolysis_rate);
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
    void SetMeteoTransport(T wind_direction,
			   T wind_speed, 
			   T pblh,
			   T ust,
			   T lmo,
			   T pressure,
			   T temperature,
			   T rain,
                           T specific_humidity);
    T GetPressure() const;
    T GetTemperature() const;
    T GetWindSpeed() const;
    T GetRain() const;
    T GetLiquidWaterContent() const;
    void SetpH(T ph);

    virtual void SetAttenuation(T attenuation);
    virtual T GetAttenuation() const;
    T GetpH() const;
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
    //
    T GetWallDryDepositionFlux(int s) const;
    T GetStreetScavengingFlux(int s) const;
    void SetMeteoDeposition(T temperature, T pressure);
    void SetMeteoScavenging(T temperature, T pressure, T rain);
    void SetStreetDryDepositionVelocity(T street_dry_deposition_velocity, int s);
    T GetStreetDryDepositionVelocity(int s) const;
    void SetWallDryDepositionVelocity(T wall_dry_deposition_velocity, int s);
    void SetStreetScavengingCoefficient(T street_scavenging_coefficient, int s);
    void SetStreetDryDepositionFlux(T street_dry_deposition_flux, int s);
    T GetStreetDryDepositionFlux(int s) const;
    T GetWallDryDepositionVelocity(int s) const;
    T GetStreetScavengingCoefficient(int s) const;
    T GetDepositionRate() const;
    void SetDepositionRate(T deposition_rate);
    
    //
    T GetTransferVelocity() const;
    void SetTransferVelocity(T transfer_velocity);
    T GetSigmaW() const;
    void SetSigmaW(T sigma_w);
    bool GetStationary() const;
    void SetStationary(bool is_stationary);

    T GetRelativeHumidity() const;
    T GetSpecificHumidity() const;
    void SetRelativeHumidity(T relative_humidity);
    //aerosol------------------------------------------

    virtual T GetStreetConcentration_aer(int s, int b) const;
    virtual T GetStreetNumberConcentration(int b) const;
    virtual void SetStreetConcentration_aer(T concentration, int s, int b);
    virtual void SetStreetNumberResuspension(T number_resuspension, int b);
    virtual void SetStreetResuspensionFactor(T resuspension_factor);
    virtual T GetStreetResuspensionFactor() const;
    virtual void SetStreetWashoffFactor(T resuspension_factor, int s);
    virtual T GetStreetWashoffFactor(int s) const;
    virtual T GetStreetNumberResuspension(int b) const;
    virtual void SetStreetWetDiameter_aer(T wet_diameter_aer, int b);
    virtual T GetStreetWetDiameter_aer(int b) const;
    virtual void SetStreetNumberConcentration(T concentration, int b);
    virtual void SetBackgroundConcentration_aer(T bg_concentration, int s, int b);
    virtual void SetEmission_aer(Array<T, 2> emission_aer);
    virtual void SetEmission_aer(T emission_aer, int s, int b);
    void SetLiquidWaterContent(T liquidwatercontent);
    virtual void SetStreetDryDepositionVelocity_aer(T street_dry_deposition_velocity, int b);
    virtual void SetWallDryDepositionVelocity_aer(T wall_dry_deposition_velocity, int b);
    virtual void SetStreetSurfaceDepositedMass_aer(T street_surface_deposited_mass_aer, int s, int b);
    virtual T GetStreetSurfaceDepositedNumber(int b) const;
    virtual void SetStreetSurfaceDepositedNumber(T street_surface_deposited_number, int b);
    virtual T GetStreetSurfaceDepositedMass_aer(int s, int b) const;
    virtual void SetStreetScavengingCoefficient_aer(T street_scavenging_coefficient, int b);
    virtual T GetStreetDryDepositionVelocity_aer(int b) const;
    virtual T GetWallDryDepositionVelocity_aer(int b) const;
    virtual T GetStreetScavengingCoefficient_aer(int b) const;
    virtual void SetStreetDryDepositionFlux_aer(T street_dry_deposition_flux, int b);
    virtual void SetWallDryDepositionFlux_aer(T wall_dry_deposition_flux, int b);
    virtual T GetStreetScavengingFluxOverCanopy_aer(int b) const;
    virtual void SetStreetScavengingFlux_aer(T street_scavenging_flux, int b);
    virtual void SetStreetNumberDryDepositionFlux(T street_number_dry_deposition_flux, int b);
    virtual void SetWallNumberDryDepositionFlux(T wall_number_dry_deposition_flux, int b);
    virtual T GetStreetNumberScavengingFluxOverCanopy(int b) const;
    virtual void SetStreetNumberScavengingFlux(T street_number_scavenging_flux, int b);
    virtual T GetStreetDryDepositionFlux_aer(int b) const;
    virtual T GetWallDryDepositionFlux_aer(int b) const;
    virtual T GetStreetScavengingFlux_aer(int b) const;
    virtual T GetMassfluxFromBackground_aer(int s, int t) const;
    virtual T GetInflowRate_aer(int s, int b) const;
    virtual T GetBackgroundConcentration_aer(int s, int b) const;
    virtual void SetMassfluxFromBackground_aer(T massflux_from_bg, int s, int b);
    virtual void SetInflowRate_aer(T inflow_rate, int s, int b);
    virtual T GetMassfluxToBackground_aer(int s, int b) const;
    virtual void SetMassfluxToBackground_aer(T massflux_to_bg, int s, int b);
    virtual T GetNumberfluxFromBackground(int b) const;
    virtual T GetNumberInflowRate(int b) const;
    virtual void SetNumberfluxFromBackground(T numberflux_from_bg, int b);
    virtual void SetNumberInflowRate(T number_inflow_rate, int b);
    virtual T GetNumberfluxToBackground(int b) const;
    virtual void SetNumberfluxToBackground(T numberflux_to_bg, int b);
    virtual T GetInitialStreetConcentration_aer(int s, int b) const;
    virtual T GetEmission_aer(int s, int b) const;
    virtual void SetMassfluxRoofToBackground_aer(T massflux_roof_to_bg_aer, int s, int b);
    virtual void SetStreetQuantityDelta_aer(T street_quantity_delta, int s, int b);
    virtual T GetInitialStreetNumberConcentration(int b) const;
    virtual T GetNumberEmission(int b) const;
    virtual T GetBackgroundNumberConcentration(int b) const;
    virtual void SetNumberfluxRoofToBackground(T massflux_roof_to_bg_aer, int b);
    virtual void SetStreetNumberQuantityDelta(T street_number_quantity_delta, int b);
    virtual void SetNumberEmission(Array<T, 1> number_emission);
    virtual void SetBackgroundNumberConcentration(T bg_number_concentration, int b);
    virtual void SetInitialStreetConcentration_aer(T concentration_aer, int s, int b);
    virtual void SetInitialStreetNumberConcentration(T number_concentration, int b);
    virtual void SetStreetRoadTraffic(T traffic_2R,
				      T traffic_HDV,
				      T traffic_PC,
				      T traffic_LCV);
    virtual T GetStreetRoadTraffic_2R() const;
    virtual T GetStreetRoadTraffic_HDV() const;
    virtual T GetStreetRoadTraffic_PC() const;
    virtual T GetStreetRoadTraffic_LCV() const;
    T GetStreetTypo() const;
    //-------------------------------------------------
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
    void SetIntersectionUST(T ustar);
    T GetLMO() const;
    T GetWindSpeed() const;
    T GetWindDirection() const;
    T GetSigmaTheta() const;
    void SetSigmaTheta(T sigma_theta);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETTRANSPORT_HXX
#endif
