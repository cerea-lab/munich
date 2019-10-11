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
  Street<T>::Street(int street_id, int begin_inter, int end_inter,
                    T length, T width, T height,
                    int ns_local, int nr_photolysis):
    street_id_(street_id), begin_inter_(begin_inter), end_inter_(end_inter),
    length_(length), width_(width), height_(height), 
    ns_local_(ns_local)
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
    photolysis_rate_.resize(nr_photolysis);
    photolysis_rate_ = 0.0;
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

  //! Returns the intersection id which the street length.
  /*!
    \return The intersection id.
  */
  template<class T>
  inline int Street<T>::GetBeginIntersectionID() const
  {
    return begin_inter_;
  }

  //! Returns the street length.
  /*!
    \return The street length (m).
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
  template<class T>
  inline void Street<T>::SetCoordinate(T longitude, T latitude)
  {
    longitude_ = longitude;
    latitude_ = latitude;
  }

  //! Returns the latitude of the street center.
  template<class T>
  inline T Street<T>::GetPhotolysisRate(int r) const
  {
    return photolysis_rate_(r);
  }

  //! Sets the street coordinate.
  template<class T>
  inline void Street<T>::SetPhotolysisRate(Array<T, 1> photolysis_rate)
  {
    photolysis_rate_ = photolysis_rate;
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


  //! Sets the meteo data.
  /*!
    \param wind_direction the wind direction.
    \param wind_speed the wind speed.
    \param pblh the PBL height.
    \param ust the friction velocity.
    \param lmo the Monin-Obukhov length
  */
  template<class T>
  inline void Street<T>::SetMeteo(T wind_direction, T wind_speed, 
                                  T pblh, T ust, T lmo)
  {
    wind_direction_ = wind_direction;
    wind_speed_ = wind_speed;
    pblh_ = pblh;
    ust_ = ust;
    lmo_ = lmo;
  }

  //! Returns the attenuation
  /*!
    \return The attenuation.
  */
  template<class T>
  inline T Street<T>::GetAttenuation() const
  {
    return attenuation_;
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
  inline void Street<T>::SetMeteoChemistry(T attenuation, T specific_humidity, 
                                           T pressure, T temperature)
  {
    attenuation_ = attenuation;
    specific_humidity_ = specific_humidity;
    pressure_ = pressure;
    temperature_ = temperature;
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

