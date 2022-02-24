#ifndef POLYPHEMUS_FILE_MODELS_STREETTRANSPORT_CXX

#include "StreetTransport.hxx"

namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  //! Main constructor.
  /*!

   */
  template<class T>
  Street<T>::Street(int street_id,
		    int begin_inter,
		    int end_inter,
		    T length,
		    T width,
		    T height,
		    int typo,
		    int ns_local):
    street_id_(street_id), begin_inter_(begin_inter), end_inter_(end_inter),
    length_(length), width_(width), height_(height), typo_(typo), ns_local_(ns_local)
  {
    emission_.resize(ns_local_);
    emission_ = 0.0;
    inflow_rate_.resize(ns_local_);
    inflow_rate_ = 0.0;
    massflux_roof_to_background_.resize(ns_local_);
    massflux_roof_to_background_ = 0.0;
    massflux_background_to_roof_.resize(ns_local_);
    massflux_background_to_roof_ = 0.0;
    massflux_to_background_.resize(ns_local_);
    massflux_to_background_ = 0.0;
    massflux_from_background_.resize(ns_local_);
    massflux_from_background_ = 0.0;
    background_concentration_.resize(ns_local_);
    background_concentration_ = 0.0;
    concentration_.resize(ns_local_);
    concentration_ = 0.0;
    initial_concentration_.resize(ns_local_);
    initial_concentration_ = 0.0;
    street_quantity_delta_.resize(ns_local_);
    street_quantity_delta_ = 0.0;
    deposition_rate_ = 0.0;
    is_stationary_ = false;
    is_low_wind_ = false;
    wind_direction_ = 0.0;
    wind_speed_ = 0.0;
    pblh_ = 0.0;
    ust_ = 0.0;
    lmo_ = 0.0;
    pressure_ = 0.0;
    temperature_ = 0.0;
    rain_ = 0.0;
    liquidwatercontent_ = 0.0;
    cloud_height_ = 0.0;
    relative_humidity_ = 0.0;
    specific_humidity_ = 0.0;
    pH_ = 4.5;
    sH_ = 1.0;

    /*** Dry deposition ***/
    street_dry_deposition_velocity_.resize(ns_local_);
    street_dry_deposition_velocity_ = 0.0;
    street_dry_deposition_flux_.resize(ns_local_);
    street_dry_deposition_flux_ = 0.0;
    wall_dry_deposition_velocity_.resize(ns_local_);
    wall_dry_deposition_velocity_ = 0.0;
    wall_dry_deposition_flux_.resize(ns_local_);
    wall_dry_deposition_flux_ = 0.0;
    roof_dry_deposition_velocity_.resize(ns_local_);
    roof_dry_deposition_velocity_ = 0.0;

    /*** Scavenging ***/
    street_scavenging_flux_overcanopy_.resize(ns_local_);
    street_scavenging_flux_overcanopy_ = 0.0;
    street_scavenging_coefficient_.resize(ns_local_);
    street_scavenging_coefficient_ = 0.0;
    street_scavenging_flux_.resize(ns_local_);
    street_scavenging_flux_ = 0.0;
  }

  //! Destructor
  template<class T>
  Street<T>::~Street()
  {
  }

  //////////////////////////////////////////
  // ACCESS METHODS FOR STREET ATTRIBUTES //
  //////////////////////////////////////////

  //! Returns the street id.
  /*!
    \return The street id.
  */
  template<class T>
  inline int Street<T>::GetStreetID() const
  {
    return street_id_;
  }

  //! Returns the intersection id where street begin.
  /*!
    \return The intersection id.
  */
  template<class T>
  inline int Street<T>::GetBeginIntersectionID() const
  {
    return begin_inter_;
  }

  //! Returns the intersection id where street end.
  /*!
    \return The intersection id.
  */
  template<class T>
  inline int Street<T>::GetEndIntersectionID() const
  {
    return end_inter_;
  }

  //! Returns the longitude of the street center.
  template<class T>
  inline T Street<T>::GetLongitude() const
  {
    return longitude_;
  }

  //! Returns the latitude of the street center.
  template<class T>
  inline T Street<T>::GetLatitude() const
  {
    return latitude_;
  }

  //! Sets the street coordinate.
  //! Returns the latitude of the street center.
  template<class T>
  inline void Street<T>::SetCoordinate(T longitude, T latitude)
  {
    longitude_ = longitude;
    latitude_ = latitude;
  }
  
  template<class T>
  inline T Street<T>::GetPhotolysisRate(int r) const
  {
    return 1;
  }

  //! Sets the photolysis rate.
  template<class T>
  inline void Street<T>::SetPhotolysisRate(Array<T, 1> photolysis_rate)
  {
    throw string("\"Street<T>::SetPhotolysisRate(Array<T, 1> photolysis_rate)\"")
      + " is not defined.";
  }

  //! Returns the street length.
  /*!
    \return The street length (m).
  */
  template<class T>
  inline T Street<T>::GetLength() const
  {
    return length_;
  }

  //! Returns the street width.
  /*!
    \return The street width (m).
  */
  template<class T>
  inline T Street<T>::GetWidth() const
  {
    return width_;
  }

  //! Returns the building height.
  /*!
    \return The building height in the street (m).
  */
  template<class T>
  inline T Street<T>::GetHeight() const
  {
    return height_;
  }

  //! Returns the building height.
  /*!
    \return The building height in the street (m).
  */
  // template<class T>
  // inline int Street<T>::GetIsTunnel() const
  // {
  //   return is_tunnel_;
  // }
  
  //! Returns the street angle.
  /*!
    \return The street angle from the begin intersection to the end intersection (rad).
  */
  template<class T>
  inline T Street<T>::GetStreetAngle() const
  {
    return street_angle_;
  }

  //! Returns the street angle.
  /*!
    \return The street angle in the intersection (rad).
  */
  template<class T>
  inline T Street<T>::GetStreetAngleIntersection() const
  {
    return street_angle_inter_;
  }

  //! Sets the street angle.
  /*!
    \param street_angle the street angle.
  */
  template<class T>
  inline void Street<T>::SetStreetAngle(T street_angle)
  {
    street_angle_ = street_angle;
  }

  //! Sets the street angle in the intersection.
  /*!
    \param street_angle the street angle.
  */
  template<class T>
  inline void Street<T>::SetStreetAngleIntersection(T street_angle)
  {
    street_angle_inter_ = street_angle;
  }

  //! Returns the wind speed in the street.
  /*!
    \return The wind speed in the street (m/s).
  */
  template<class T>
  inline T Street<T>::GetStreetWindSpeed() const
  {
    return ustreet_;
  }

  //! Sets the wind speed in the street.
  /*!
    \param ustreet the wind speed in the street.
  */
  template<class T>
  inline void Street<T>::SetStreetWindSpeed(T ustreet)
  {
    ustreet_ = ustreet;
  }

  //! Returns the friction velocity.
  /*!
    \return The friction velocity in the street (m/s).
  */
  template<class T>
  inline T Street<T>::GetStreetUstar() const
  {
    return ust_;
  }

  //! Sets the friction velocity.
  /*!
    \param ustar_street the friction velocity.
  */
  template<class T>
  inline void Street<T>::SetStreetUstar(T ustar_street)
  {
    ust_ = ustar_street;
  }

  //! Returns the canopy volume.
  /*!
    \return The canopy volume (m3).
  */
  template<class T>
  inline T Street<T>::GetStreetVolume() const
  {
    return width_ * length_ * height_;
  }

  //! Returns the mass quantity.
  /*!
    \return The mass quantity in the street (ug).
  */
  template<class T>
  inline T Street<T>::GetStreetQuantity(int s) const
  {
    T volume = width_ * length_ * height_;
    return concentration_(s) * volume;
  }

  //! Returns the concentration.
  /*!
    \return The concentration in the street (ug/m3).
  */
  template<class T>
  inline T Street<T>::GetStreetConcentration(int s) const
  {
    return concentration_(s);
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetConcentration(T concentration, int s)
  {
    concentration_(s) = concentration;
  }

  //! Returns the initial concentration.
  /*!
    \return The initial concentration in the street (ug/m3).
  */
  template<class T>
  inline T Street<T>::GetInitialStreetConcentration(int s) const
  {
    return initial_concentration_(s);
  }

  //! Sets the initial concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetInitialStreetConcentration(T concentration, int s)
  {
    initial_concentration_(s) = concentration;
  }

  //! Returns the mass quantity change.
  /*!
    \return The mass quantity changen in the street (ug).
  */
  template<class T>
  inline T Street<T>::GetStreetQuantityDelta(int s) const
  {
    return street_quantity_delta_(s);
  }


  //! Sets the mass quantity change.
  /*!
    \param concentration the mass quantity changeconcentration.
  */
  template<class T>
  inline void Street<T>::SetStreetQuantityDelta(T street_quantity_delta, int s)
  {
    street_quantity_delta_(s) = street_quantity_delta;
  }

  //! Returns the background concentration.
  /*!
    \return The background concentration for the street (ug/m3).
  */
  template<class T>
  inline T Street<T>::GetBackgroundConcentration(int s) const
  {
    return background_concentration_(s);
  }

  //! Sets the background concentration.
  /*!
    \param background_concentration the background concentration.
  */
  template<class T>
  inline void Street<T>::SetBackgroundConcentration(T background_concentration, int s)
  {
    background_concentration_(s) = background_concentration;
  }

  //! Returns the emission rate.
  /*!
    \return The emission rate in the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetEmission(int s) const
  {
    return emission_(s);
  }

  //! Sets the emission rate.
  /*!
    \param emission The emission rate in the street (ug/s).
  */
  template<class T>
  inline void Street<T>::SetEmission(Array<T, 1> emission)
  {
    emission_ = emission;
  }

    //! Returns the emission rate.
  /*!
    \return The emission rate in the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetWindSpeed() const
  {
    return wind_speed_;
  }

  //! Returns the wind direction.
  /*!
    \return The wind direction in the street (radian).
  */
  template<class T>
  inline T Street<T>::GetWindDirection() const
  {
    return wind_direction_;
  }

  //! Returns the planetary boundary layer height (PBL height).
  /*!
    \return The PBL height (m).
  */
  template<class T>
  inline T Street<T>::GetPBLH() const
  {
    return pblh_;
  }

  //! Returns the Monin-Obukhov length.
  /*!
    \return The Monin-Obukhov length in the street (m).
  */
  template<class T>
  inline T Street<T>::GetLMO() const
  {
    return lmo_;
  }

  //! Returns rain in the street.
  /*!
    \return rain in the street (mm).
  */
  template<class T>
  inline T Street<T>::GetRain() const
  {
    return rain_;
  }


  //! Sets the meteo data.
  /*!
    \param wind_direction the wind direction.
    \param wind_speed the wind speed.
    \param pblh the PBL height.
    \param ust the friction velocity.
    \param lmo the Monin-Obukhov length
    \param pressure the pressure
    \param temperature the temperature
  */
  template<class T>
  inline void Street<T>::SetMeteoTransport(T wind_direction,
					   T wind_speed, 
					   T pblh,
					   T ust,
					   T lmo,
                                           T pressure,
					   T temperature,
					   T rain,
					   T specific_humidity)
  {
    wind_direction_ = wind_direction;
    wind_speed_ = wind_speed;
    pblh_ = pblh;
    ust_ = ust;
    lmo_ = lmo;
    pressure_ = pressure;
    temperature_ = temperature;
    rain_ = rain;
    specific_humidity_ = specific_humidity;
  }

  //! Returns the attenuation
  /*!
    \return The attenuation.
  */
  template<class T>
  inline T Street<T>::GetAttenuation() const
  {
    return 1;
  }

  //! Returns the specific humidity
  /*!
    \return The specific humidity.
  */
  template<class T>
  inline T Street<T>::GetSpecificHumidity() const
  {
    return specific_humidity_;
  }

  //! Returns the pressure
  /*!
    \return The pressure.
  */
  template<class T>
  inline T Street<T>::GetPressure() const
  {
    return pressure_;
  }

  //! Returns the temperature
  /*!
    \return The temperature.
  */
  template<class T>
  inline T Street<T>::GetTemperature() const
  {
    return temperature_;
  }

  //! Sets the meteo data.
  /*!
    \param attenuation the attenuation.
    \param specific_humidity the specific humidity.
    \param pressure the pressure.
    \param temperature the temperature.
  */
  template<class T>
  inline void Street<T>::SetAttenuation(T attenuation)
  {
    throw string("\"Street<T>::SetAttenuation(T attenuation)\"")
      + " is not defined.";
  }

  //! 
  /*!
    \return (1/s)
  */
  template<class T>
  inline T Street<T>::GetStreetScavengingCoefficient(int s) const
  {
    return street_scavenging_coefficient_(s);
  }

  //! 
  /*!
    \param 1/s
  */
  template<class T>
  inline void Street<T>::SetStreetScavengingCoefficient(T street_scavenging_coefficient, int s)
  {
    street_scavenging_coefficient_(s) = street_scavenging_coefficient;
  }  

  //! Returns the inflow rate.
  /*!
    \return The inflow rate to the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetInflowRate(int s) const
  {
    return inflow_rate_(s);
  }

  //! Sets the inflow rate.
  /*!
    \param inflow_rate the inflow rate (ug/s).
  */
  template<class T>
  inline void Street<T>::SetInflowRate(T inflow_rate, int s)
  {
    inflow_rate_(s) = inflow_rate;
  }

  //! Returns the outgoing volume flux.
  /*!
    \return The outgoing volume flux (m3/s).
  */
  template<class T>
  inline T Street<T>::GetOutgoingFlux() const
  {
    return outgoing_flux_;
  }

  //! Sets the outgoing volume flux.
  /*!
    \param outgoing_flux the outgoing volume flux (m3/s).
  */
  template<class T>
  inline void Street<T>::SetOutgoingFlux(T outgoing_flux)
  {
    outgoing_flux_ = outgoing_flux;
  }

  //! Returns the incoming volume flux from other streets.
  /*!
    \return The incoming volume flux (m3/s).
  */
  template<class T>
  inline T Street<T>::GetIncomingFlux() const
  {
    return incoming_flux_;
  }

  //! Sets the incoming volume flux.
  /*!
    \param incoming_flux the outgoing volume flux (m3/s).
  */
  template<class T>
  inline void Street<T>::SetIncomingFlux(T incoming_flux)
  {
    incoming_flux_ = incoming_flux;
  }


  //! Returns the mass flux to the atmosphere.
  /*!
    \return The mass flux to the atmosphere (ug/s).
  */
  template<class T>
  inline T Street<T>::GetMassfluxToBackground(int s) const
  {
    return massflux_to_background_(s);
  }

  //! Sets th outflow rate to the atmosphere.
  /*!
    \param massflux_to_background the outflow rate (ug/s).
  */
  template<class T>
  inline void Street<T>::SetMassfluxToBackground(T massflux_to_background, int s)
  {
    massflux_to_background_(s) = massflux_to_background;
  }

  //! Returns the mass flux from the atmosphere.
  /*!
    \return The mass flux from the atmosphere (ug/s).
  */
  template<class T>
  inline T Street<T>::GetMassfluxFromBackground(int s) const
  {
    return massflux_from_background_(s);
  }

  //! Sets the mass flux from the atmosphere.
  /*!
    \param massflux_from_background the mass flux from the atmosphere (ug/s).
  */
  template<class T>
  inline void Street<T>::SetMassfluxFromBackground(T massflux_from_background, int s)
  {
    massflux_from_background_(s) = massflux_from_background;
  }

  //! Returns the mass flux to background through the roof.
  /*!
    \return The mass flux to background through the roof (ug/s).
  */
  template<class T>
  inline T Street<T>::GetMassfluxRoofToBackground(int s) const
  {
    return massflux_roof_to_background_(s);
  }

  //! Sets  the mass flux to background through the roof.
  /*!
    \param massflux_roof_to_background the mass flux to background through the roof (ug/s).
  */
  template<class T>
  inline void Street<T>::SetMassfluxRoofToBackground(T massflux_roof_to_background, int s)
  {
    massflux_roof_to_background_(s) = massflux_roof_to_background;
  }

  //! Returns the mass flux from the atmosphere through the roof.
  /*!
    \return The mass flux from the atmosphere through the roof (ug/s).
  */
  template<class T>
  inline T Street<T>::GetMassfluxBackgroundToRoof(int s) const
  {
    return massflux_background_to_roof_(s);
  }

  //! Sets the mass flux from the atmosphere through the roof.
  /*!
    \param massflux_background_to_roof the mass flux from the atmosphere through the roof (ug/s).
  */
  template<class T>
  inline void Street<T>::SetMassfluxBackgroundToRoof(T massflux_background_to_roof, int s)
  {
    massflux_background_to_roof_(s) = massflux_background_to_roof;
  }

  //! Returns the deposition rate.
  /*!
    \return The deposition rate (ug/s).
  */
  template<class T>
  inline T Street<T>::GetDepositionRate() const
  {
    return deposition_rate_;
  }

  //! Sets the deposition rate.
  /*!
    \param deposition_rate the deposition rate (ug/s).
  */
  template<class T>
  inline void Street<T>::SetDepositionRate(T deposition_rate)
  {
    deposition_rate_ = deposition_rate;
  }

  
  template<class T>
  inline T Street<T>::GetStreetDryDepositionVelocity(int s) const
  {
    return street_dry_deposition_velocity_(s);
  }

    template<class T>
  inline void Street<T>::SetStreetDryDepositionVelocity(T street_dry_deposition_velocity,
                                                        int s)
  {
    street_dry_deposition_velocity_(s) = street_dry_deposition_velocity;
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetDryDepositionFlux(T street_dry_deposition_flux, int s)
  {
    street_dry_deposition_flux_(s) = street_dry_deposition_flux;
  }

  
  //! 
  /*!
    \return (m3/s)
  */
  template<class T>
  inline T Street<T>::GetWallDryDepositionVelocity(int s) const
  {
    return wall_dry_deposition_velocity_(s);
  }

    //! 
  /*!
    \param 
  */
  template<class T>
  inline void Street<T>::SetWallDryDepositionVelocity(T wall_dry_deposition_velocity,
                                                      int s)
  {
    wall_dry_deposition_velocity_(s) = wall_dry_deposition_velocity;
  }

  //! Returns the transfer velocity at roof.
  /*!
    \return The transfer velocity at roof (m/s).
  */
  template<class T>
  inline T Street<T>::GetTransferVelocity() const
  {
    return transfer_velocity_;
  }

  //! Sets the transfer velocity at roof. 
  /*!
    \param transfer_velocity the transfer velocity at roof (m/s).
  */
  template<class T>
  inline void Street<T>::SetTransferVelocity(T transfer_velocity)
  {
    transfer_velocity_ = transfer_velocity;
  }

  //! Returns the standard deviation of the vertical velocity.
  /*!
    \return The standard deviation of the vertical velocity (m/s).
  */
  template<class T>
  inline T Street<T>::GetSigmaW() const
  {
    return sigma_w_;
  }

  //! Sets the standard deviation of the vertical velocity.
  /*!
    \param sigma_w the standard deviation of the vertical velocity (m/s).
  */
  template<class T>
  inline void Street<T>::SetSigmaW(T sigma_w)
  {
    sigma_w_ = sigma_w;
  }

  //! Returns the stationarity.
  /*!
    \return The stationarity.
  */
  template<class T>
  inline bool Street<T>::GetStationary() const
  {
    return is_stationary_;
  }

  //! Sets the stationarity.
  /*!
    \param is_stationary the stationarity.
  */
  template<class T>
  inline void Street<T>::SetStationary(bool is_stationary)
  {
    is_stationary_ = is_stationary;
  }

  
  //aerosol ----------------------------------------------
  template<class T>
  inline T Street<T>::GetRelativeHumidity() const
  {
    return relative_humidity_;
  }

  //! Sets the street coordinate.
  template<class T>
  inline void Street<T>::SetRelativeHumidity(T relative_humidity)
  {
    relative_humidity_ = relative_humidity;
  }

  //! Returns the aerosol concentration.
  /*!
    \return The concentration in the street (ug/m3).
  */
  template<class T>
  inline T Street<T>::GetStreetConcentration_aer(int s, int b) const
  {
    return 1;
  }

  //! Returns the number concentration.
  /*!
    \return The number concentration in the street (#/m3).
  */
  template<class T>
  inline T Street<T>::GetStreetNumberConcentration(int b) const
  {
    return 1;
  }

  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetConcentration_aer(T concentration, int s, int b)
  {
    throw string("\"StreetAerosol<T>::SetStreetConcentration_aer(T concentration, int s, int b)\"")
      + " is not defined.";
  }

  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetRoadTraffic(T traffic_2R,
						     T traffic_HDV,
						     T traffic_PC,
						     T traffic_LCV)
  {
    throw string("\"StreetAerosol<T>::SetStreetRoadTraffic(T traffic_2R, T traffic_HDV, T traffic_PC, T traffic_LCV)\"")
      + " is not defined.";
  }

  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline T Street<T>::GetStreetRoadTraffic_2R() const
  {
    throw string("\"GetStreetRoadTraffic_2R()\"")
      + " is not defined.";
  }

  template<class T>
  inline T Street<T>::GetStreetRoadTraffic_HDV() const
  {
    throw string("\"GetStreetRoadTraffic_HDV()\"")
      + " is not defined.";
  }

  template<class T>
  inline T Street<T>::GetStreetRoadTraffic_PC() const
  {
    throw string("\"GetStreetRoadTraffic_PC()\"")
      + " is not defined.";
  }

  template<class T>
  inline T Street<T>::GetStreetRoadTraffic_LCV() const
  {
    throw string("\"GetStreetRoadTraffic_LCV()\"")
      + " is not defined.";
  }

    //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetNumberResuspension(T number_resuspension, int b)
  {
    throw string("\"StreetAerosol<T>::SetStreetNumberResuspension(T number_resuspension, int b)\"")
      + " is not defined.";
  }

  
  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetResuspensionFactor(T resuspension_factor)
  {
    throw string("\"StreetAerosol<T>::SetStreetResuspensionFactor(T resuspension_factor)\"")
      + " is not defined.";
  }

  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline T Street<T>::GetStreetResuspensionFactor() const
  {
    return 1;
  }

  
  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetWashoffFactor(T washoff_factor, int s)
  {
    throw string("\"StreetAerosol<T>::SetStreetWashoffFactor(T washoff_factor, int s)\"")
      + " is not defined.";
  }

  
  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline T Street<T>::GetStreetWashoffFactor(int s) const
  {
    return 1;
  }

  
  //! Sets the aerosol concentration.
  /*!
    \param concentration_aer the aerosol concentration.
  */
  template<class T>
  inline T Street<T>::GetStreetNumberResuspension(int b) const
  {
    return 1;
  }

  
  //! Sets the aerosol wet diameter.
  template<class T>
  inline void Street<T>::SetStreetWetDiameter_aer(T wet_diameter_aer, int b)
  {
    throw string("\"StreetAerosol<T>::SetStreetWetDiameter_aer(T wet_diameter_aer, int b)\"")
      + " is not defined.";
  }

  
  //! Sets the number concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetNumberConcentration(T concentration, int b)
  {
    throw string("\"Street<T>::SetStreetNumberConcentration(T concentration, int s, int b)\"")
      + " is not defined.";
  }

  
  //! Sets the background concentration.
  /*!
    \param bg_concentration the background concentration.
  */
  template<class T>
  inline void Street<T>::SetBackgroundConcentration_aer(T bg_concentration, int s, int b)
  {
    throw string("\"Street<T>::SetBackgroundConcentration_aer(int, int)\"")
      + " is not defined.";
  }

  
  //! Sets the aerosol emission rate.
  /*!
    \param emission_aer The aerosol emission rate in the street (ug/s).
  */
  template<class T>
  inline void Street<T>::SetEmission_aer(Array<T, 2> emission_aer)
  {
    throw string("\"Street<T>::SetEmission_aer(Array<T, 2> emission_aer)\"")
      + " is not defined.";
  }
  
  //! Sets the aerosol emission rate.
  /*!
    \param emission_aer The aerosol emission rate in the street (ug/s).
  */
  template<class T>
  inline void Street<T>::SetEmission_aer(T emission_aer, int s, int b)
  {
    throw string("\"Street<T>::SetEmission_aer(T emission_aer, int s, int b)\"")
      + " is not defined.";
  }

  template<class T>
  inline void Street<T>::SetLiquidWaterContent(T liquidwatercontent)
  {
    liquidwatercontent_ = liquidwatercontent;
  }

  //! Returns the emission rate.
  /*!
    \return The emission rate in the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetLiquidWaterContent() const
  {
    return liquidwatercontent_;
  }

  //! Returns the emission rate.
  /*!
    \return The emission rate in the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetStreetTypo() const
  {
    return typo_;
  }


    //! 
  /*!
    \param 
  */
  template<class T>
  inline void Street<T>::SetStreetDryDepositionVelocity_aer(T street_dry_deposition_velocity, int b)
  {
    throw string("\"Street<T>::SetStreetDryDepositionVelocity_aer(T street_dry_deposition_velocity, int b)\"")
      + " is not defined.";
  }

  
  //! 
  /*!
    \param 
  */
  template<class T>
  inline void Street<T>::SetWallDryDepositionVelocity_aer(T wall_dry_deposition_velocity, int b)
  {
    throw string("\"Street<T>::SetWallDryDepositionVelocity_aer(T wall_dry_deposition_velocity, int b)\"")
      + " is not defined.";
  }

  
  //! 
  /*!
    \param 
  */
  template<class T>
  inline void Street<T>::SetStreetSurfaceDepositedMass_aer(T street_surface_deposited_mass, int s, int b)
  {
    throw string("\"Street<T>::SetStreetSurfaceDepositedMass_aer(T street_surface_deposited_mass_aer, int s, int b)\"")
      + " is not defined.";
  }

  
  //! 
  /*!
    \param 
  */
  template<class T>
  inline void Street<T>::SetStreetSurfaceDepositedNumber(T street_surface_deposited_number, int b)
  {
    throw string("\"Street<T>::SetStreetSurfaceDepositedMass_aer(T street_surface_deposited_number, int b)\"")
      + " is not defined.";
  }

  template<class T>
  inline T Street<T>::GetStreetSurfaceDepositedNumber(int b) const
  {
    return 1;
  }

  template<class T>
  inline T Street<T>::GetStreetSurfaceDepositedMass_aer(int s, int b) const
  {
    return 1;
  }

  
  //! 
  /*!
    \param 1/s
  */
  template<class T>
  inline void Street<T>::SetStreetScavengingCoefficient_aer(T street_scavenging_coefficient, int b)
  {
    throw string("\"Street<T>::SetStreetScavengingCoefficient_aer(T street_scavenging_coefficient, int b)\"")
      + " is not defined.";
  }

  template<class T>
  inline T Street<T>::GetStreetDryDepositionVelocity_aer(int b) const
  {
    return 1;
  }

  //! 
  /*!
    \return (m3/s)
  */
  template<class T>
  inline T Street<T>::GetWallDryDepositionVelocity_aer(int b) const
  {
    return 1;
  }

  
  //! 
  /*!
    \return (1/s)
  */
  template<class T>
  inline T Street<T>::GetStreetScavengingCoefficient_aer(int b) const
  {
    return 1;
  }

  
  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetDryDepositionFlux_aer(T street_dry_deposition_flux, int b)
  {
    throw string("\"Street<T>::SetStreetDryDepositionFlux_aer(T street_dry_deposition_flux, int b)\"")
      + " is not defined.";
  }

  
  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetWallDryDepositionFlux_aer(T wall_dry_deposition_flux, int b)
  {
    throw string("\"Street<T>::SetWallDryDepositionFlux_aer(T wall_dry_deposition_flux, int b)\"")
      + " is not defined.";
  }

  
  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T Street<T>::GetStreetScavengingFluxOverCanopy_aer(int b) const
  {
    return 1;
  }

  
  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetScavengingFlux_aer(T street_scavenging_flux, int b)
  {
    throw string("\"Street<T>::SetStreetScavengingFlux_aer(T street_scavenging_flux, int b)\"")
      + " is not defined.";
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetNumberDryDepositionFlux(T street_number_dry_deposition_flux, int b)
  {
    throw string("\"Street<T>SetStreetNumberDryDepositionFlux(T street_number_dry_deposition_flux, int b)\"")
      + " is not defined.";
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetWallNumberDryDepositionFlux(T wall_number_dry_deposition_flux, int b)
  {
    throw string("\"Street<T>SetWallNumberDryDepositionFlux(T wall_number_dry_deposition_flux, int b)\"")
      + " is not defined.";
  }

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T Street<T>::GetStreetNumberScavengingFluxOverCanopy(int b) const
  {
    throw string("\"Street<T>::GetStreetNumberScavengingFluxOverCanopy(int b)\"")
      + " is not defined.";
  }

  //! Sets the concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetStreetNumberScavengingFlux(T street_number_scavenging_flux, int b)
  {
    throw string("\"Street<T>::SetStreetNumberScavengingFlux(T street_number_scavenging_flux, int b)\"")
      + " is not defined.";
  }

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T Street<T>::GetStreetDryDepositionFlux_aer(int b) const
  {
    return 1;
  }

  
  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T Street<T>::GetWallDryDepositionFlux_aer(int b) const
  {
    return 1;
  }

  //! Returns the street dry deposition flux.
  /*!
    \return The street dry deposition flux in the street (ug/m2/s).
  */
  template<class T>
  inline T Street<T>::GetStreetScavengingFlux_aer(int b) const
  {
    return 1;
  }

  
  //! Returns the background concentration.
  /*!
    \return The background concentration for the street (ug/m3).
  */
  template<class T>
  inline T Street<T>::GetBackgroundConcentration_aer(int s, int b) const
  {
    return 1;
  }

  //! Sets the mass flux from the atmosphere.
  /*!
    \param massflux_from_bg the mass flux from the atmosphere (ug/s).
  */
  template<class T>
  inline void Street<T>::SetMassfluxFromBackground_aer(T massflux_from_bg, int s, int b)
  {
    throw string("\"Street<T>::SetMassfluxFromBackground_aer(T massflux_from_bg, int s, int b)\"")
      + " is not defined.";
  }  
  
  //! Sets the inflow rate.
  /*!
    \param inflow_rate the inflow rate (ug/s).
  */
  template<class T>
  inline void Street<T>::SetInflowRate_aer(T inflow_rate, int s, int b)
  {
    throw string("\"Street<T>::SetInflowRate_aer(T inflow_rate, int s, int b)\"")
      + " is not defined.";
  }  
  
  //! Returns the mass flux to the atmosphere.
  /*!
    \return The mass flux to the atmosphere (ug/s).
  */
  template<class T>
  inline T Street<T>::GetMassfluxToBackground_aer(int s, int b) const
  {
    return 1;
  }  
  
  //! Sets th outflow rate to the atmosphere.
  /*!
    \param massflux_to_bg the outflow rate (ug/s).
  */
  template<class T>
  inline void Street<T>::SetMassfluxToBackground_aer(T massflux_to_bg, int s, int b)
  {
    throw string("\"Street<T>::SetMassfluxToBackground_aer(T massflux_to_bg, int s, int b)\"")
      + " is not defined.";
  }  
  
  // //Number
  
  //! Returns the mass flux from the atmosphere.
  /*!
    \return The mass flux from the atmosphere (#/s).
  */
  template<class T>
  inline T Street<T>::GetNumberfluxFromBackground(int b) const
  {
    return 1;
  } 
  
  
  //! Returns the inflow rate.
  /*!
    \return The inflow rate to the street (#/s).
  */
  template<class T>
  inline T Street<T>::GetNumberInflowRate(int b) const
  {
    return 1;
  }

  //! Sets the mass flux from the atmosphere.
  /*!
    \param massflux_from_bg the mass flux from the atmosphere (#/s).
  */
  template<class T>
  inline void Street<T>::SetNumberfluxFromBackground(T numberflux_from_bg, int b)
  {
    throw string("\"Street<T>::SetNumberfluxFromBackground(T numberflux_from_bg, int b)\"")
      + " is not defined.";
  }  
  
  //! Sets the inflow rate.
  /*!
    \param inflow_rate the inflow rate (#/s).
  */
  template<class T>
  inline void Street<T>::SetNumberInflowRate(T number_inflow_rate, int b)
  {
    throw string("\"Street<T>::SetNumberInflowRate(T number_inflow_rate, int b)\"")
      + " is not defined.";
  }  
  
  //! Returns the mass flux to the atmosphere.
  /*!
    \return The mass flux to the atmosphere (#/s).
  */
  template<class T>
  inline T Street<T>::GetNumberfluxToBackground(int b) const
  {
    return 1;
  }  
  
  //! Sets th outflow rate to the atmosphere.
  /*!
    \param massflux_to_bg the outflow rate (#/s).
  */
  template<class T>
  inline void Street<T>::SetNumberfluxToBackground(T numberflux_to_bg, int b)
  {
    throw string("\"Street<T>::SetNumberInflowRate(T number_inflow_rate, int b)\"")
      + " is not defined.";    
  }

  //Concentration - mass

  //! Returns the initial concentration.
  /*!
    \return The initial concentration in the street (ug/m3).
  */
  template<class T>
  inline T Street<T>::GetInitialStreetConcentration_aer(int s, int b) const
  {
    return 1;
  }
  
  //! Returns the emission rate.
  /*!
    \return The emission rate in the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetEmission_aer(int s, int b) const
  {
    return 1;
  }  
  
  //! Sets  the mass flux to background through the roof.
  /*!
    \param massflux_roof_to_bg the mass flux to background through the roof (ug/s).
  */
  template<class T>
  inline void Street<T>::SetMassfluxRoofToBackground_aer(T massflux_roof_to_bg_aer, int s, int b)
  {
    throw string("\"Street<T>::SetMassfluxRoofToBackground_aer(T massflux_roof_to_bg_aer, int s, int b)\"")
      + " is not defined.";    
  }  
  
  //! Sets the mass quantity change.
  /*!
    \param concentration the mass quantity changeconcentration.
  */
  template<class T>
  inline void Street<T>::SetStreetQuantityDelta_aer(T street_quantity_delta, int s, int b)
  {
     throw string("\"Street<T>::SetStreetQuantityDelta_aer(T street_quantity_delta, int s, int b)\"")
      + " is not defined.";  
  }  
  
  //! Returns the initial concentration.
  /*!
    \return The initial concentration in the street (#/m3).
  */
  template<class T>
  inline T Street<T>::GetInitialStreetNumberConcentration(int b) const
  {
    return 1;
  }  
  
  //! Returns the emission rate.
  /*!
    \return The emission rate in the street (#/s).
  */
  template<class T>
  inline T Street<T>::GetNumberEmission(int b) const
  {
    return 1;
  }  
  
  //! Returns the background concentration.
  /*!
    \return The background concentration for the street (ug/m3).
  */
  template<class T>
  inline T Street<T>::GetBackgroundNumberConcentration(int b) const
  {
    return 1;
  }  
  
  //! Sets  the mass flux to background through the roof.
  /*!
    \param massflux_roof_to_bg the mass flux to background through the roof (ug/s).
  */
  template<class T>
  inline void Street<T>::SetNumberfluxRoofToBackground(T massflux_roof_to_bg_aer, int b)
  {
    throw string("\"Street<T>::SetNumberfluxRoofToBackground(T massflux_roof_to_bg_aer, int b)\"")
      + " is not defined.";    
  }
  //! Sets the mass quantity change.
  /*!
    \param concentration the mass quantity changeconcentration.
  */
  template<class T>
  inline void Street<T>::SetStreetNumberQuantityDelta(T street_number_quantity_delta, int b)
  {
     throw string("\"Street<T>::SetStreetNumberQuantityDelta(T street_number_quantity_delta, int b)\"")
      + " is not defined.";  
  }  
  
  //! Sets the aerosol emission rate.
  /*!
    \param emission_aer The aerosol emission rate in the street (ug/s).
  */
  template<class T>
  inline void Street<T>::SetNumberEmission(Array<T, 1> number_emission)
  {
    throw string("\"Street<T>::SetNumberEmission(Array<T, 1> number_emission)\"")
      + " is not defined.";
  }  
  
  //! Sets the background concentration.
  /*!
    \param bg_concentration the background concentration.
  */
  template<class T>
  inline void Street<T>::SetBackgroundNumberConcentration(T bg_number_concentration, int b)
  {
    throw string("\"Street<T>::SetBackgroundNumberConcentration(T bg_number_concentration, int b)\"")
      + " is not defined.";    
  }  
  
  //! Sets the initial concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetInitialStreetConcentration_aer(T concentration_aer, int s, int b)
  {
    throw string("\"Street<T>::SetInitialStreetConcentration_aer(T concentration_aer, int s, int b)\"")
      + " is not defined.";   
  }  
  
  //! Sets the initial concentration.
  /*!
    \param concentration the concentration.
  */
  template<class T>
  inline void Street<T>::SetInitialStreetNumberConcentration(T number_concentration, int b)
  {
    throw string("\"Street<T>::SetInitialStreetNumberConcentration(T number_concentration, int b)\"")
      + " is not defined."; 
  }

  //! Sets the aerosol wet diameter.
  template<class T>
  inline T Street<T>::GetStreetWetDiameter_aer(int b) const
  {
    return 1;
  }  

  
  //! Returns the mass flux from the atmosphere.
  /*!
//     \return The mass flux from the atmosphere (ug/s).
//   */
  template<class T>
  inline T Street<T>::GetMassfluxFromBackground_aer(int s, int t) const
  {
    return 1;
  }

  
  //! Returns the inflow rate.
  /*!
    \return The inflow rate to the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetInflowRate_aer(int s, int b) const
  {
    return 1;
  }

  //! Sets the meteo data to calculate dry deposition.
  /*!
    \param temperature is the ambient temperature.
    \param pressure is the ambient pressure.
  */
  template<class T>
  inline void Street<T>::SetpH(T ph)
  {
    pH_ = ph;
  }

  //! Returns the emission rate.
  /*!
    \return The emission rate in the street (ug/s).
  */
  template<class T>
  inline T Street<T>::GetpH() const
  {
    return pH_;
  }

  //! Sets the mixing length factor sH (-)
  template<class T>
  inline void Street<T>::SetsH(T sH)
  {
    sH_ = sH;
  }

  //! Returns the mixing length factor sH
  template<class T>
  inline T Street<T>::GetsH() const
  {
    return sH_;
  }
 
  //--------------------------------------------------------


  ////////////////////////
  // Intersection class //
  ////////////////////////

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  //! Main constructor.
  /*!

   */
  template<class T>
  Intersection<T>::Intersection(int id, T x, T y, int nstreet, 
                                Array<int, 1> street_list, bool is_virtual,
                                Array<T, 2> flux_matrix, Array<T, 2> gaussian_matrix):
    id_(id), x_(x), y_(y), nstreet_(nstreet), street_list_(street_list), 
    is_virtual_(is_virtual), flux_matrix_(flux_matrix), 
    gaussian_matrix_(gaussian_matrix)
  {
    wind_direction_ = 0.0;
    wind_speed_ = 0.0;
    pblh_ = 0.0;
    ust_ = 0.0;
    lmo_ = 0.0;
    sigma_theta_ = 0.0;
  }

  //! Destructor
  template<class T>
  Intersection<T>::~Intersection()
  {
  }

  ////////////////////////////////////////////////
  // ACCESS METHODS FOR INTERSECTION ATTRIBUTES //
  ////////////////////////////////////////////////

  //! Returns the intersection id.
  /*!
    \return The intersection id.
  */
  template<class T>
  inline int Intersection<T>::GetID() const
  {
    return id_;
  }

  //! Returns the intersection longitude.
  /*!
    \return The intersection longitude.
  */
  template<class T>
  inline T Intersection<T>::GetX() const
  {
    return x_;
  }

  //! Returns the intersection latitude.
  /*!
    \return The intersection latitude.
  */
  template<class T>
  inline T Intersection<T>::GetY() const
  {
    return y_;
  }

  //! Returns the number of streets which are connected to the intersection.
  /*!
    \return The number of streets.
  */
  template<class T>
  inline T Intersection<T>::GetNStreet() const
  {
    return nstreet_;
  }

  //! Returns the list of street id which are connected to the intersection.
  /*!
    \return The list of street id.
  */
  template<class T>
  inline Array<int, 1>& Intersection<T>::GetStreetList() 
  {
    return street_list_;
  }

  //! Is the intersection a end node of the street network?
  /*!
    \return .
  */
  template<class T>
  inline bool Intersection<T>::IsVirtual() 
  {
    return is_virtual_;
  }


  //! Returns the flux matrix.
  /*!
    \return The flux matrix.
  */
  template<class T>
  inline Array<T, 2>& Intersection<T>::GetFluxMatrix() 
  {
    return flux_matrix_;
  }

  //! Sets the flux matrix.
  /*!
    \param flux_matrix the flux matrix.
  */
  template<class T>
  inline void Intersection<T>::SetFluxMatrix(Array<T, 2> flux_matrix)
  {
    flux_matrix_ = flux_matrix;
  }

  //! Returns the gaussian matrix.
  /*!
    \return The gaussian matrix.
  */
  template<class T>
  inline Array<T, 2>& Intersection<T>::GetGaussianMatrix() 
  {
    return gaussian_matrix_;
  }

  //! Sets the gaussian matrix.
  /*!
    \param gaussian matrix the gaussian matrix.
  */
  template<class T>
  inline void Intersection<T>::SetGaussianMatrix(Array<T, 2> gaussian_matrix)
  {
    gaussian_matrix_ = gaussian_matrix;
  }

  //! Sets the meteo data for the intersection.
  /*!
    \param wind_direction the wind direction.
    \param wind_speed the wind speed.
    \param pblh the PBL height.
    \param ust the friction velocity.
    \param lmo the Monin-Obukhov length.
  */
  template<class T>
  inline void Intersection<T>::SetMeteo(T wind_direction, T wind_speed, 
                                        T pblh, T ust, T lmo)
  {
    wind_direction_ = wind_direction;
    wind_speed_ = wind_speed;
    pblh_ = pblh;
    ust_ = ust;
    lmo_ = lmo;
  }

  //! Returns the PBL height.
  /*!
    \return The PBL height (m).
  */
  template<class T>
  inline T Intersection<T>::GetPBLH() const
  {
    return pblh_;
  }

  //! Returns the friction velocity.
  /*!
    \return The friction velocity (m/s).
  */
  template<class T>
  inline T Intersection<T>::GetUST() const
  {
    return ust_;
  }

  template<class T>
  inline void Intersection<T>::SetIntersectionUST(T ustar)
  {
    ust_ = ustar;
  }


  //! Returns the Monin-Obukhov length.
  /*!
    \return The Monin-Obukhov length (m).
  */
  template<class T>
  inline T Intersection<T>::GetLMO() const
  {
    return lmo_;
  }

  //! Returns the wind speed.
  /*!
    \return The wind speed (m/s).
  */
  template<class T>
  inline T Intersection<T>::GetWindSpeed() const
  {
    return wind_speed_;
  }

  //! Returns the wind direction.
  /*!
    \return The wind direction (radian).
  */
  template<class T>
  inline T Intersection<T>::GetWindDirection() const
  {
    return wind_direction_;
  }

  //! Returns the standard deviation of the horizontal wind direction.
  /*!
    \return The standard deviation of the horizontal wind direction.
  */
  template<class T>
  inline T Intersection<T>::GetSigmaTheta() const
  {
    return sigma_theta_;
  }

  //! Sets the standard deviation of the horizontal wind direction.
  /*!
    \param sigma_theta the standard deviation of the horizontal wind direction.
  */
  template<class T>
  inline void Intersection<T>::SetSigmaTheta(T sigma_theta)
  {
    sigma_theta_ = sigma_theta;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETTRANSPORT_CXX
#endif

