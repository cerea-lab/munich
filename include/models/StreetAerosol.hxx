#ifndef POLYPHEMUS_FILE_MODELS_STREETAEROSOL_HXX

#include "StreetChemistry.cxx"

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
  class StreetAerosol: public StreetChemistry<T>
  {

  protected:

    //! Number of aerosol species.
    int ns_local_aer_;    
    
    //! Number of aerosol size sections.
    int nsize_local_;   

    //! Aerosol emission rate (ug/s). 
    Array<T, 2> emission_aer_; 

    //! Aerosol inflow rate to the street (ug/s).
    Array<T, 2> inflow_rate_aer_;
    
    //! Outflow aerosol rate to the atmosphere (ug/s)
    Array<T, 2> massflux_roof_to_bg_aer_;    

    //! List of aerosol species concentrations (ug/m3).
    Array<T, 2> concentration_aer_; 
    
    //! List of aerosol species number concentrations (#/m3).
    Array<T, 1> number_concentration_;    

    //! List of background aerosol species concentrations (ug/m3).
    Array<T, 2> bg_concentration_aer_;
    
    //! List of aerosol wet diameters (m).
    Array<T, 1> wet_diameter_aer_;
    
    //! Inflow rate from the atmosphere to the canyon through the intersection (ug/s).
    Array<T, 2> massflux_from_bg_aer_;  
    
    //! Outflow rate to the atmosphere through the intersection (ug/s)
    Array<T, 2> massflux_to_bg_aer_;
    
    //*** Resuspension
    T traffic_2R_;
    T traffic_HDV_;
    T traffic_PC_;
    T traffic_LCV_;

    //! Particles resuspension (mug/s)
    T resuspension_factor_;
    Array<T, 1> washoff_factor_;
    T pm10_resuspension_;

    //! Particle number resuspension (#/s)
    Array<T, 1> number_resuspension_;
    
    //! Particles tyre emission (mug/s)
    Array<T, 2> tyre_emission_;

    //! Particles brake emission (mug/s)
    Array<T, 2> brake_emission_;

    //! Particles road emission (mug/s)
    Array<T, 2> road_emission_; 
    
    //Number
    
    //! Inflow rate from the atmosphere to the canyon through the intersection (#/s).
    Array<T, 1> numberflux_from_bg_;   
    
    //! Inflow rate to the street (#/s).
    Array<T, 1> number_inflow_rate_;    
    
    //! List of background species concentrations (#/m3).
    Array<T, 1> bg_number_concentration_;    
    
    //! Outflow rate to the atmosphere through the intersection (#/s)
    Array<T, 1> numberflux_to_bg_; 
    
    //Concentration - mass
    
    //! List of species initial concentrations (ug/m3).
    Array<T, 2> initial_concentration_aer_;    
    
    //! List of species mass change (ug).
    Array<T, 2> street_quantity_delta_aer_;    
    
    //Concentration - number
   
    //! List of number initial concentrations (#/m3).
    Array<T, 1> initial_number_concentration_;    

    //! Aerosol emission rate (#/s). 
    Array<T, 1> number_emission_;     
    
    //! Outflow aerosol rate to the atmosphere (#/s)
    Array<T, 1> numberflux_roof_to_bg_;   
    
    //! List of number change (#).
    Array<T, 1> street_number_quantity_delta_;   

    //! Scavenging
    //! .
    Array<T, 1> street_scavenging_flux_aer_;   
    //! Outflow rate to the atmosphere through the intersection (ug/s)
    Array<T, 1> street_scavenging_coefficient_aer_;  
    //!
    Array<T, 1> street_scavenging_flux_overcanopy_aer_;   
    Array<T, 1> street_number_scavenging_flux_overcanopy_;
    Array<T, 1> street_number_scavenging_flux_;

    //! Deposition
    Array<T, 2> street_surface_deposited_mass_aer_;
    Array<T, 1> street_surface_deposited_number_;
    Array<T, 1> street_dry_deposition_flux_aer_;
    Array<T, 1> wall_dry_deposition_flux_aer_;
    Array<T, 1> tree_dry_deposition_flux_aer_;
    Array<T, 1> street_dry_deposition_velocity_aer_;
    Array<T, 1> wall_dry_deposition_velocity_aer_;
    Array<T, 1> tree_dry_deposition_velocity_aer_;
    Array<T, 1> street_number_dry_deposition_flux_;
    Array<T, 1> wall_number_dry_deposition_flux_;
    Array<T, 1> tree_number_dry_deposition_flux_;

    T relative_humidity_;
    T street_surface_water_;

  public:

     /*** Constructor and destructor ***/

    StreetAerosol(int street_id,
		  int begin_inter,
		  int end_inter, 
		  T length,
		  T width,
		  T height,
		  int typo,
		  int ns_local,
		  int nr_photolysis,
		  int ns_local_aer,
		  int nsize_local);
    virtual ~StreetAerosol();   

    /*** Methods ***/

    T GetStreetConcentration_aer(int s, int b) const;
    T GetRelativeHumidity() const;
    void SetRelativeHumidity(T relative_humidity);
    T GetStreetNumberConcentration(int b) const;
    T GetStreetWetDiameter_aer(int b) const;
    void SetStreetWetDiameter_aer(T wet_diameter_aer, int b);
    void SetStreetConcentration_aer(T concentration_aer, int s, int b);
    void SetStreetNumberConcentration(T number_concentration, int b);
    T GetInitialStreetConcentration_aer(int s, int b) const;
    void SetStreetQuantityDelta_aer(T street_quantity_delta, int s, int b);
    void SetBackgroundConcentration_aer(T bg_concentration_aer, int s, int b);
    void SetEmission_aer(Array<T, 2> emission_aer);
    void SetEmission_aer(T emission_aer, int s, int b);
    T GetEmission_aer(int s, int b) const;
    //Changes start here
    T GetBackgroundConcentration_aer(int s, int b) const;    
    T GetInflowRate_aer(int s, int b) const;    
    T GetMassfluxFromBackground_aer(int s, int b) const;
    void SetMassfluxFromBackground_aer(T massflux_from_bg, int s, int b);
    void SetInflowRate_aer(T inflow_rate, int s, int b);
    T GetMassfluxToBackground_aer(int s, int b) const;
    void SetMassfluxToBackground_aer(T massflux_to_bg, int s, int b);
    
    //Number
    T GetNumberfluxFromBackground(int b) const;
    T GetNumberInflowRate(int b) const;
    T GetBackgroundNumberConcentration(int b) const;
    void SetNumberfluxFromBackground(T numberflux_from_bg, int b);
    void SetNumberInflowRate(T number_inflow_rate, int b);
    T GetNumberfluxToBackground(int b) const;
    void SetNumberfluxToBackground(T numberflux_to_bg, int b);
    
    //Mass concentration
    void SetMassfluxRoofToBackground_aer(T massflux_roof_to_bg_aer, int s, int b);
    //Number concentration
    T GetInitialStreetNumberConcentration(int b) const;
    T GetNumberEmission(int b) const;
    
    void SetNumberfluxRoofToBackground(T numberflux_roof_to_bg, int b);
    void SetStreetNumberQuantityDelta(T street_number_quantity_delta, int b);
    
    void SetNumberEmission(Array<T, 1> number_emission);
    void SetBackgroundNumberConcentration(T bg_number_concentration, int b);
    
    void SetInitialStreetNumberConcentration(T number_concentration, int b);
    void SetInitialStreetConcentration_aer(T concentration_aer, int s, int b); 
    
    Data<T, 2>& GetStreetConcentration_aer();
    Data<T, 1>& GetStreetNumberConcentration();
    
    //Scavenging
    T GetStreetDryDepositionVelocity_aer(int b) const;
    T GetStreetScavengingFlux_aer(int b) const;
    void SetStreetScavengingCoefficient_aer(T street_scavenging_coefficient_aer, int b);
    T GetStreetScavengingCoefficient_aer(int b) const;
    T GetStreetScavengingFluxOverCanopy_aer(int b) const;
    void SetStreetScavengingFlux_aer(T street_scavenging_flux, int b);
    T GetStreetNumberScavengingFluxOverCanopy(int b) const;
    void SetStreetNumberScavengingFlux(T street_number_scavenging_flux, int b);

    T GetStreetSurfaceWater() const;
    void SetStreetSurfaceWater(T street_surface_water);
    //Deposition
    T GetStreetSurfaceDepositedMass_aer(int s, int b) const;
    void SetStreetSurfaceDepositedMass_aer(T street_surface_deposited_mass_aer, int s, int b);
    T GetStreetSurfaceDepositedNumber(int b) const;
    void SetStreetSurfaceDepositedNumber(T street_surface_deposited_number, int b);
    T GetStreetDryDepositionFlux_aer(int b) const;
    T GetWallDryDepositionFlux_aer(int b) const;
    T GetTreeDryDepositionFlux_aer(int b) const;
    void SetStreetDryDepositionVelocity_aer(T street_dry_deposition_velocity, int b);
    void SetWallDryDepositionVelocity_aer(T wall_dry_deposition_velocity, int b);
    T GetWallDryDepositionVelocity_aer(int b) const;
    void SetTreeDryDepositionVelocity_aer(T tree_dry_deposition_velocity, int b);
    T GetTreeDryDepositionVelocity_aer(int b) const;
    void SetStreetDryDepositionFlux_aer(T street_dry_deposition_flux, int b);
    void SetWallDryDepositionFlux_aer(T wall_dry_deposition_flux, int b);
    void SetTreeDryDepositionFlux_aer(T street_dry_deposition_flux, int b);
    void SetStreetNumberDryDepositionFlux(T street_number_dry_deposition_flux, int b);
    void SetWallNumberDryDepositionFlux(T wall_number_dry_deposition_flux, int b);
    void SetTreeNumberDryDepositionFlux(T tree_number_dry_deposition_flux, int b);
    void SetStreetWashoffFactor(T washoff_factor, int s);
    T GetStreetWashoffFactor(int s) const;
    
    //Resuspension
    void SetStreetResuspensionFactor(T resuspension_factor);
    T GetStreetResuspensionFactor() const;
    void SetStreetPM10Resuspension(T pm10_resuspension);
    T GetStreetPM10Resuspension() const;
    void SetStreetNumberResuspension(T number_resuspension, int b);
    T GetStreetNumberResuspension(int b) const;
    void SetStreetTyreEmission(T tyre_emission, int s, int b);
    T GetStreetTyreEmission(int s, int b) const;
    void SetStreetBrakeEmission(T brake_emission, int s, int b);
    T GetStreetBrakeEmission(int s, int b) const;
    void SetStreetRoadEmission(T road_emission, int s, int b);
    T GetStreetRoadEmission(int s, int b) const;

    //Road traffic
    void SetStreetRoadTraffic(T traffic_2R,
			      T traffic_HDV,
			      T traffic_PC,
			      T traffic_LCV);
    
    T GetStreetRoadTraffic_2R() const;
    T GetStreetRoadTraffic_HDV() const;
    T GetStreetRoadTraffic_PC() const;
    T GetStreetRoadTraffic_LCV() const;
    
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETAEROSOL_HXX
#endif
