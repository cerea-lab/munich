#ifndef POLYPHEMUS_FILE_MODELS_STREETAEROSOL_CXX

#include "StreetAerosol.hxx"

namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  //! Main constructor.
  /*!

   */
  template<class T>
  StreetAerosol<T>::StreetAerosol(int street_id,
				  int begin_inter,
				  int end_inter,
				  T length,
				  T width,
				  T height,
				  int typo,
				  int ns_local,
				  int nr_photolysis,
				  int ns_local_aer,
				  int nsize_local):
    StreetChemistry<T> (street_id,
			begin_inter,
			end_inter,
			length,
			width,
			height,
			typo,
			ns_local,
			nr_photolysis),
    ns_local_aer_(ns_local_aer),
    nsize_local_(nsize_local)
  {
    emission_aer_.resize(ns_local_aer_, nsize_local_);
    emission_aer_ = 0.0;
    inflow_rate_aer_.resize(ns_local_aer_, nsize_local_);
    inflow_rate_aer_ = 0.0;    
    massflux_roof_to_bg_aer_.resize(ns_local_aer_, nsize_local_);
    massflux_roof_to_bg_aer_ = 0.0;    
    bg_concentration_aer_.resize(ns_local_aer_, nsize_local_);
    bg_concentration_aer_ = 0.0;
    concentration_aer_.resize(ns_local_aer_, nsize_local_);
    concentration_aer_ = 0.0;
    number_concentration_.resize(nsize_local_);
    number_concentration_ = 0.0;
    massflux_from_bg_aer_.resize(ns_local_aer_, nsize_local_);
    massflux_from_bg_aer_ = 0.0;
    massflux_to_bg_aer_.resize(ns_local_aer_, nsize_local_);
    massflux_to_bg_aer_ = 0.0;

    numberflux_from_bg_.resize(nsize_local_);
    numberflux_from_bg_ = 0.0;
    number_inflow_rate_.resize(nsize_local_);
    number_inflow_rate_ = 0.0;
    bg_number_concentration_.resize(nsize_local_);
    bg_number_concentration_ = 0.0;
    numberflux_to_bg_.resize(nsize_local_);
    numberflux_to_bg_ = 0.0; 

    initial_concentration_aer_.resize(ns_local_aer_, nsize_local_);
    initial_concentration_aer_ = 0.0; 
    street_quantity_delta_aer_.resize(ns_local_aer_, nsize_local_);
    street_quantity_delta_aer_ = 0.0; 

    initial_number_concentration_.resize(nsize_local_);
    initial_number_concentration_ = 0.0;
    number_emission_.resize(nsize_local_);
    number_emission_ = 0.0;
    numberflux_roof_to_bg_.resize(nsize_local_);
    numberflux_roof_to_bg_ = 0.0;
    street_number_quantity_delta_.resize(nsize_local_);
    street_number_quantity_delta_ = 0.0;

    street_scavenging_coefficient_aer_.resize(nsize_local_);
    street_scavenging_coefficient_aer_ = 0.0;
    street_scavenging_flux_overcanopy_aer_.resize(nsize_local_);
    street_scavenging_flux_overcanopy_aer_ = 0.0;

    street_surface_deposited_mass_aer_.resize(ns_local_aer_, nsize_local_);
    street_surface_deposited_mass_aer_ = 0.0;
    street_surface_deposited_number_.resize(nsize_local_);
    street_surface_deposited_number_ = 0.0;

    street_surface_water_ = 0.0; 

    street_dry_deposition_velocity_aer_.resize(nsize_local_);
    street_dry_deposition_velocity_aer_ = 0.0;
    wall_dry_deposition_velocity_aer_.resize(nsize_local_);
    wall_dry_deposition_velocity_aer_ = 0.0;
    street_dry_deposition_flux_aer_.resize(nsize_local_);
    street_dry_deposition_flux_aer_ = 0.0;
    wall_dry_deposition_flux_aer_.resize(nsize_local_);
    wall_dry_deposition_flux_aer_ = 0.0;
    street_number_dry_deposition_flux_.resize(nsize_local_);
    street_number_dry_deposition_flux_ = 0.0;
    wall_number_dry_deposition_flux_.resize(nsize_local_);
    wall_number_dry_deposition_flux_ = 0.0;
    street_number_scavenging_flux_overcanopy_.resize(nsize_local_);
    street_number_scavenging_flux_overcanopy_ = 0.0;
    street_number_scavenging_flux_.resize(nsize_local_);
    street_number_scavenging_flux_ = 0.0;
    street_scavenging_flux_aer_.resize(nsize_local_);
    street_scavenging_flux_aer_ = 0.0;

    wet_diameter_aer_.resize(nsize_local_);
    wet_diameter_aer_ = 0.0;

    washoff_factor_.resize(ns_local_aer_);
    washoff_factor_ = 0.0;

    //Nh_week_ = 168;
    // nb_hdv_hour_.resize(Nh_week_);
    // nb_hdv_hour_ = 0.0;         
    // nb_ldv_hour_.resize(Nh_week_);
    // nb_ldv_hour_ = 0.0;         
    // u_hdv_.resize(Nh_week_);
    // u_hdv_ = 0.0;         
    // u_ldv_.resize(Nh_week_);
    // u_ldv_ = 0.0;
    traffic_2R_ = 0.0;
    traffic_HDV_ = 0.0;
    traffic_PC_ = 0.0;
    traffic_LDV_ = 0.0;

    resuspension_factor_ = 0.0;
    pm10_resuspension_ = 0.0;
    number_resuspension_.resize(nsize_local_);
    number_resuspension_ = 0.0;
    tyre_emission_.resize(ns_local_aer_, nsize_local_);
    tyre_emission_ = 0.0;
    brake_emission_.resize(ns_local_aer_, nsize_local_);
    brake_emission_ = 0.0;
    road_emission_.resize(ns_local_aer_, nsize_local_);
    road_emission_ = 0.0;

    relative_humidity_ = 0.0;
  }

  //! Destructor
  template<class T>
  StreetAerosol<T>::~StreetAerosol()
  {
  }

  //////////////////////////////////////////
  // ACCESS METHODS FOR STREET ATTRIBUTES //
  //////////////////////////////////////////
  
  /*!
    \return The aerosol concentration in the street (ug/m3).
  */  
  template<class T>
  inline T StreetAerosol<T>::GetStreetConcentration_aer(int s, int b) const
  {
    return concentration_aer_(s,b);
  } 
  
  //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  Data<T, 2>& StreetAerosol<T>::GetStreetConcentration_aer()
  {
    return concentration_aer_;
  }  

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetNumberScavengingFluxOverCanopy(int b) const
  {
    return street_number_scavenging_flux_overcanopy_(b);
  }

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetRelativeHumidity() const
  {
    return relative_humidity_;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetRelativeHumidity(T relative_humidity)
  {
    relative_humidity_ = relative_humidity;
  }
  
  /*!
    \return The aerosol number concentration in the street (ug/m3).
  */  
  template<class T>
  inline T StreetAerosol<T>::GetStreetNumberConcentration(int b) const
  {
    return number_concentration_(b);
  }  
  
  //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  Data<T, 1>& StreetAerosol<T>::GetStreetNumberConcentration()
  {
    return number_concentration_;
  }    

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetWetDiameter_aer(T wet_diameter_aer, int b)
  {
    wet_diameter_aer_(b) = wet_diameter_aer;
  }
  
  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetResuspensionFactor(T resuspension_factor)
  {
    resuspension_factor_ = resuspension_factor;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetRoadTraffic(T traffic_2R,
						     T traffic_HDV,
						     T traffic_PC,
						     T traffic_LDV)
  {
    traffic_2R_ = traffic_2R;
    traffic_HDV_ = traffic_HDV;
    traffic_PC_ = traffic_PC;
    traffic_LDV_ = traffic_LDV;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetWashoffFactor(T washoff_factor, int s)
  {
    washoff_factor_(s) = washoff_factor;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetPM10Resuspension(T pm10_resuspension)
  {
    pm10_resuspension_ = pm10_resuspension;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetNumberResuspension(T number_resuspension, int b)
  {
    number_resuspension_(b) = number_resuspension;
  }

    /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetTyreEmission(T tyre_emission, int s, int b)
  {
    tyre_emission_(s,b) = tyre_emission;
  }

    //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetTyreEmission(int s, int b) const
  {
    return tyre_emission_(s, b);
  }

      /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetBrakeEmission(T brake_emission, int s, int b)
  {
    brake_emission_(s,b) = brake_emission;
  }

    //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetBrakeEmission(int s, int b) const
  {
    return brake_emission_(s, b);
  }

      /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetRoadEmission(T road_emission, int s, int b)
  {
    road_emission_(s,b) = road_emission;
  }

    //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetRoadEmission(int s, int b) const
  {
    return road_emission_(s, b);
  } 

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetResuspensionFactor() const
  {
    return resuspension_factor_;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetRoadTraffic_2R() const
  {
    return traffic_2R_;
  }
  
  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetRoadTraffic_HDV() const
  {
    return traffic_HDV_;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetRoadTraffic_PC() const
  {
    return traffic_PC_;
  }

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetRoadTraffic_LDV() const
  {
    return traffic_LDV_;
  }
  
  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetWashoffFactor(int s) const
  {
    return washoff_factor_(s);
  }    

    /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetPM10Resuspension() const
  {
    return pm10_resuspension_;
  }  
    /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetNumberResuspension(int b) const
  {
    return number_resuspension_(b);
  }    

  /*!
    \param wet_diameter_aer the aerosol wet diameter.
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetWetDiameter_aer(int b) const
  {
    return wet_diameter_aer_(b);
  }  
  
  /*!
    \param aerosol_concentration the aerosol concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetConcentration_aer(T concentration_aer, int s, int b)
  {
    concentration_aer_(s,b) = concentration_aer;
  }
  
    template<class T>
  inline void StreetAerosol<T>::SetStreetNumberConcentration(T number_concentration, int b)
  {
    number_concentration_(b) = number_concentration;
  }

  //! Sets the aerosol background concentration.
  /*!
    \param bg_concentration_aer the aerosol background concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetBackgroundConcentration_aer(T bg_concentration_aer, int s, int b)
  {
    bg_concentration_aer_(s, b) = bg_concentration_aer;
  }
  
  //! Returns the aerosol emission rate.
  /*!
    \return The aerosol emission rate in the street (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetEmission_aer(int s, int b) const
  {
    return emission_aer_(s, b);
  }  
  
    //! Sets the aerosol emission rate.
  /*!
    \param emission The aerosol emission rate in the street (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetEmission_aer(Array<T, 2> emission_aer)
  {
    emission_aer_ = emission_aer;
  }

  //! Sets the aerosol emission rate.
  /*!
    \param emission The aerosol emission rate in the street (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetEmission_aer(T emission_aer, int s, int b)
  {
    emission_aer_(s, b) = emission_aer;
  }
 
  //! 
  /*!
    \param 1/s
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetScavengingCoefficient_aer(T street_scavenging_coefficient_aer, int b)
  {
    street_scavenging_coefficient_aer_(b) = street_scavenging_coefficient_aer;
  }

  //! Returns the aerosol background concentration.
  /*!
    \return The aerosol background concentration for the street (ug/m3).
  */
  template<class T>
  inline T StreetAerosol<T>::GetBackgroundConcentration_aer(int s, int b) const
  {
    return bg_concentration_aer_(s, b);
  }  

  //! Returns the aerosol inflow rate.
  /*!
    \return The aerosol inflow rate to the street (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetInflowRate_aer(int s, int b) const
  {
    return inflow_rate_aer_(s, b);
  }

  //! Returns the mass flux from the atmosphere.
  /*!
    \return The mass flux from the atmosphere (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetMassfluxFromBackground_aer(int s, int b) const
  {
    return massflux_from_bg_aer_(s, b);
  }
  
  //! Sets the mass flux from the atmosphere.
  /*!
    \param massflux_from_bg the mass flux from the atmosphere (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetMassfluxFromBackground_aer(T massflux_from_bg, int s, int b)
  {
    massflux_from_bg_aer_(s, b) = massflux_from_bg;
  }  
  
  //! Sets the inflow rate.
  /*!
    \param inflow_rate the inflow rate (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetInflowRate_aer(T inflow_rate, int s, int b)
  {
    inflow_rate_aer_(s, b) = inflow_rate;
  }  
  
  //! Returns the mass flux to the atmosphere.
  /*!
    \return The mass flux to the atmosphere (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetMassfluxToBackground_aer(int s, int b) const
  {
    return massflux_to_bg_aer_(s, b);
  }  
  
  //! Sets th outflow rate to the atmosphere.
  /*!
    \param massflux_to_bg the outflow rate (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetMassfluxToBackground_aer(T massflux_to_bg, int s, int b)
  {
    massflux_to_bg_aer_(s, b) = massflux_to_bg;
  }  
  
  //Number
  
  //! Returns the mass flux from the atmosphere.
  /*!
    \return The mass flux from the atmosphere (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetNumberfluxFromBackground(int b) const
  {
    return numberflux_from_bg_(b);
  }  
  
    //! Returns the inflow rate.
  /*!
    \return The inflow rate to the street (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetNumberInflowRate(int b) const
  {
    return number_inflow_rate_(b);
  }  
  
  //! Returns the background concentration.
  /*!
    \return The background concentration for the street (ug/m3).
  */
  template<class T>
  inline T StreetAerosol<T>::GetBackgroundNumberConcentration(int b) const
  {
    return bg_number_concentration_(b);
  }    
  
  //! Sets the mass flux from the atmosphere.
  /*!
    \param numberflux_from_bg the mass flux from the atmosphere (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetNumberfluxFromBackground(T numberflux_from_bg, int b)
  {
    numberflux_from_bg_(b) = numberflux_from_bg;
  }  
  
  //! Sets the inflow rate.
  /*!
    \param inflow_rate the inflow rate (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetNumberInflowRate(T number_inflow_rate, int b)
  {
    number_inflow_rate_(b) = number_inflow_rate;
  }  
  
  //! Returns the mass flux to the atmosphere.
  /*!
    \return The mass flux to the atmosphere (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetNumberfluxToBackground(int b) const
  {
    return numberflux_to_bg_(b);
  }  
  
  //! Sets th outflow rate to the atmosphere.
  /*!
    \param numberflux_to_bg the outflow rate (#/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetNumberfluxToBackground(T numberflux_to_bg, int b)
  {
    numberflux_to_bg_(b) = numberflux_to_bg;
  }  
  
  //! Returns the initial concentration.
  /*!
    \return The initial concentration in the street (ug/m3).
  */
  template<class T>
  inline T StreetAerosol<T>::GetInitialStreetConcentration_aer(int s, int b) const
  {
    return initial_concentration_aer_(s, b);
  }  

    //! Sets  the mass flux to background through the roof.
  /*!
    \param massflux_roof_to_bg the mass flux to background through the roof (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetMassfluxRoofToBackground_aer(T massflux_roof_to_bg_aer, int s, int b)
  {
    massflux_roof_to_bg_aer_(s, b) = massflux_roof_to_bg_aer;
  }  
  
  //! Sets the mass quantity change.
  /*!
    \param concentration the mass quantity changeconcentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetQuantityDelta_aer(T street_quantity_delta_aer, int s, int b)
  {
    street_quantity_delta_aer_(s, b) = street_quantity_delta_aer;
  }  
  
//Concentration Number
  //! Returns the initial concentration.
  /*!
    \return The initial concentration in the street (ug/m3).
  */
  template<class T>
  inline T StreetAerosol<T>::GetInitialStreetNumberConcentration(int b) const
  {
    return initial_number_concentration_(b);
  }
  
  //! Returns the emission rate.
  /*!
    \return The emission rate in the street (ug/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetNumberEmission(int b) const
  {
    return number_emission_(b);
  }  

  template<class T>
  inline T StreetAerosol<T>::GetStreetSurfaceWater() const
  {
    return street_surface_water_;
  }  

  //    */
  template<class T>
  inline void StreetAerosol<T>::SetStreetSurfaceWater(T street_surface_water)
  {
    street_surface_water_ = street_surface_water;
  }  
  
  template<class T>
  inline T StreetAerosol<T>::GetStreetSurfaceDepositedMass_aer(int s, int b) const
  {
    return street_surface_deposited_mass_aer_(s, b);
  }  

  //    */
  template<class T>
  inline void StreetAerosol<T>::SetStreetSurfaceDepositedMass_aer(T street_surface_deposited_mass_aer, int s, int b)
  {
    street_surface_deposited_mass_aer_(s, b) = street_surface_deposited_mass_aer;
  }  


  template<class T>
  inline T StreetAerosol<T>::GetStreetSurfaceDepositedNumber(int b) const
  {
    return street_surface_deposited_number_(b);
  }  

  //    */
  template<class T>
  inline void StreetAerosol<T>::SetStreetSurfaceDepositedNumber(T street_surface_deposited_number, int b)
  {
    street_surface_deposited_number_(b) = street_surface_deposited_number;
  }  
  
  template<class T>
  inline T StreetAerosol<T>::GetStreetDryDepositionVelocity_aer(int b) const
  {
    return street_dry_deposition_velocity_aer_(b);
  }  
  
  //! Sets  the mass flux to background through the roof.
  /*!
    \param numberflux_roof_to_bg the number flux to background through the roof (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetNumberfluxRoofToBackground(T numberflux_roof_to_bg, int b)
  {
    numberflux_roof_to_bg_(b) = numberflux_roof_to_bg;
  }  
  
  //! Sets the mass quantity change.
  /*!
    \param concentration the mass quantity changeconcentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetNumberQuantityDelta(T street_number_quantity_delta, int b)
  {
    street_number_quantity_delta_(b) = street_number_quantity_delta;
  }   
  
  
  //! Sets the emission rate.
  /*!
    \param emission The emission rate in the street (ug/s).
  */
  template<class T>
  inline void StreetAerosol<T>::SetNumberEmission(Array<T, 1> number_emission)
  {
    number_emission_ = number_emission;
  } 
  
  //! Sets the background concentration.
  /*!
    \param bg_concentration the background concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetBackgroundNumberConcentration(T bg_number_concentration, int b)
  {
    bg_number_concentration_(b) = bg_number_concentration;
  }  
  
  //! Sets the initial concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetInitialStreetConcentration_aer(T concentration_aer, int s, int b)
  {
    initial_concentration_aer_(s, b) = concentration_aer;
  }  
  
  //! Sets the initial concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetInitialStreetNumberConcentration(T number_concentration, int b)
  {
    initial_number_concentration_(b) = number_concentration;
  }    

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetScavengingFlux_aer(int b) const
  {
    return street_scavenging_flux_aer_(b);
  }

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetScavengingFluxOverCanopy_aer(int b) const
  {
    return street_scavenging_flux_overcanopy_aer_(b);
  }


  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetNumberScavengingFlux(T street_number_scavenging_flux, int b)
  {
    street_number_scavenging_flux_(b) = street_number_scavenging_flux;
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetScavengingFlux_aer(T street_scavenging_flux, int b)
  {
    street_scavenging_flux_aer_(b) = street_scavenging_flux;
  }

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetDryDepositionFlux_aer(int b) const
  {
    return street_dry_deposition_flux_aer_(b);
  }

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T StreetAerosol<T>::GetWallDryDepositionFlux_aer(int b) const
  {
    return wall_dry_deposition_flux_aer_(b);
  }

  //! 
  /*!
    \return (m3/s)
  */
  template<class T>
  inline T StreetAerosol<T>::GetWallDryDepositionVelocity_aer(int b) const
  {
    return wall_dry_deposition_velocity_aer_(b);
  }


  //! 
  /*!
    \return (1/s)
  */
  template<class T>
  inline T StreetAerosol<T>::GetStreetScavengingCoefficient_aer(int b) const
  {
    return street_scavenging_coefficient_aer_(b);
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetDryDepositionFlux_aer(T street_dry_deposition_flux, int b)
  {
    street_dry_deposition_flux_aer_(b) = street_dry_deposition_flux;
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetNumberDryDepositionFlux(T street_number_dry_deposition_flux, int b)
  {
    street_number_dry_deposition_flux_(b) = street_number_dry_deposition_flux;
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetWallNumberDryDepositionFlux(T wall_number_dry_deposition_flux, int b)
  {
    wall_number_dry_deposition_flux_(b) = wall_number_dry_deposition_flux;
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void StreetAerosol<T>::SetWallDryDepositionFlux_aer(T wall_dry_deposition_flux, int b)
  {
    wall_dry_deposition_flux_aer_(b) = wall_dry_deposition_flux;
  }

  //! 
  /*!
    \param 
  */
  template<class T>
  inline void StreetAerosol<T>::SetStreetDryDepositionVelocity_aer(T street_dry_deposition_velocity, int b)
  {
    street_dry_deposition_velocity_aer_(b) = street_dry_deposition_velocity;
  }

  //! 
  /*!
    \param 
  */
  template<class T>
  inline void StreetAerosol<T>::SetWallDryDepositionVelocity_aer(T wall_dry_deposition_velocity, int b)
  {
    wall_dry_deposition_velocity_aer_(b) = wall_dry_deposition_velocity;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETAEROSOL_CXX
#endif

