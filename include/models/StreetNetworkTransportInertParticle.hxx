#ifndef POLYPHEMUS_FILE_MODELS_STREETNETWORKTRANSPORTINERTPARTICLE_HXX

#include "BaseModel.cxx"
#include "StreetTransportInertParticle.cxx"

namespace Polyphemus
{


  using namespace std;
  using namespace blitz;

  //////////////////////
  // FORTRAN FUNCTION //
  //////////////////////


#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#ifdef POLYPHEMUS_DOUBLE_UNDERSCORE
  //#define _compute_dry_deposition compute_dry_deposition__
#define _compute_scavenging_coefficient_aer	\
  compute_scavenging_coefficient_aer__
#else
  //#define _compute_dry_deposition compute_dry_deposition_
#define _compute_scavenging_coefficient_aer	\
  compute_scavenging_coefficient_aer_
#endif

  extern "C"
  {
    void _compute_scavenging_coefficient_aer(int*, double*, double*, double*,
					     double*, double*, double*,
					     double*, double*);
  }

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

    static const T pi, earth_radius, karman, nu, Pr;

    /*** Configurations ***/

    string option_transfer, option_ustreet, option_method;
    T sub_delta_t_min;
    bool is_stationary;
    //! Output configuration.
    string output_config;
    string output_dir;
    bool text_file;

    /*** Domain ***/
    //! 2D grid for streets.
    RegularGrid<T> GridST2D;
    
    //! Grid for intersections.
    RegularGrid<T> GridINT2D;

    //! 2D grid for species.
    RegularGrid<T> GridS2D;

    /*** Emission ***/

    int Nt_emis;
    int Ns_emis;
    
    //! Emission rate at current date.
    Data<T, 2> Emission_i;    
    //! Emission rate at next date.
    Data<T, 2> Emission_f;    
    //! Emission rate buffer.
    Data<T, 2> FileEmission_i;
    //! Emission rate buffer.
    Data<T, 2> FileEmission_f;
    
    //! List of species with emissions.
    vector<string> species_list_emis;
    
    /*** Deposition ***/

    string option_wind_profile;
    //! List of species with deposition velocities.
    vector<string> species_list_dep;
    //! Number of species with deposition velocities.
    //int Ns_dep;
    //! Dry deposition fluxes at current date.
    Data<T, 2> StreetDryDepositionFlux;
    //! Dry wall deposition fluxes at current date.
    Data<T, 2> WallDryDepositionFlux;
    //! Dry deposition rate at current date.
    Data<T, 2> StreetDryDepositionRate;
    //! Dry wall deposition rate at current date.
    Data<T, 2> WallDryDepositionRate;
    
    /*** Scavenging ***/

    //! List of species with scavenging.
    vector<string> species_list_scav;
    //! Number of species with scavenging.
    //int Ns_scav;
    //! Scavenging coefficient at current date.
    Data<T, 2> ScavengingCoefficient;
    //! Scavenging fluxes at current date.
    Data<T, 2> StreetScavengingFlux;
    //! Scavenging rate at current date.
    Data<T, 2> StreetScavengingRate;
    //!  Scavenging rain threshold (mm / h).
    T scavenging_rain_threshold;
    
    /*** Street data ***/

    string file_street;
    int total_nstreet;
    Array<int, 1> id_street, begin_inter, end_inter;
    Array<T, 1> length, width, height;
    Array<int, 1> is_tunnel;
    //! Pointer to the current street.
    typename vector<Street<T>* >::iterator current_street;
    //! Street list.
    vector<Street<T>*> StreetVector, StreetVectorInter;

    /*** Background concentration data ***/

    int Nt_background;
    int Ns_background;

    //! Background concentration at current date.
    Data<T, 2> Background_i;    
    //! Background concentration at next date.
    Data<T, 2> Background_f;    
    //! Background concentration buffer.
    Data<T, 2> FileBackground_i;
    //! Background concetration buffer.
    Data<T, 2> FileBackground_f;

    //! List of species with background concentrations.
    vector<string> species_list_background;
    T cell_volume;

    /*** Meteorological data ***/

    int Nt_meteo; 
    T ustreet_min; // Minimum wind speed in the streets (m/s)
    
    //! Rain at current date.
    Data<T, 1> Rain_i;    
    //! Rain at next date.
    Data<T, 1> Rain_f;    
    //! Rain buffer.
    Data<T, 1> FileRain_i;
    //! Rain buffer.
    Data<T, 1> FileRain_f;

    //! Temperature at current date.
    Data<T, 1> Temperature_i;    
    //! Temperature at next date.
    Data<T, 1> Temperature_f;    
    //! Temperature buffer.
    Data<T, 1> FileTemperature_i;
    //! Temperature buffer.
    Data<T, 1> FileTemperature_f;

    //! Pressure at current date.
    Data<T, 1> Pressure_i;    
    //! Pressure at next date.
    Data<T, 1> Pressure_f;    
    //! Pressure buffer.
    Data<T, 1> FilePressure_i;
    //! Pressure buffer.
    Data<T, 1> FilePressure_f;

    //! WindDirection at current date.
    Data<T, 1> WindDirection_i;    
    //! WindDirection at next date.
    Data<T, 1> WindDirection_f;    
    //! WindDirection buffer.
    Data<T, 1> FileWindDirection_i;
    //! WindDirection buffer.
    Data<T, 1> FileWindDirection_f;

    //! WindSpeed at current date.
    Data<T, 1> WindSpeed_i;    
    //! WindSpeed at next date.
    Data<T, 1> WindSpeed_f;    
    //! WindSpeed buffer.
    Data<T, 1> FileWindSpeed_i;
    //! WindSpeed buffer.
    Data<T, 1> FileWindSpeed_f;

    //! PBLH at current date.
    Data<T, 1> PBLH_i;    
    //! PBLH at next date.
    Data<T, 1> PBLH_f;    
    //! PBLH buffer.
    Data<T, 1> FilePBLH_i;
    //! PBLH buffer.
    Data<T, 1> FilePBLH_f;

    //! UST at current date.
    Data<T, 1> UST_i;    
    //! UST at next date.
    Data<T, 1> UST_f;    
    //! UST buffer.
    Data<T, 1> FileUST_i;
    //! UST buffer.
    Data<T, 1> FileUST_f;

    //! LMO at current date.
    Data<T, 1> LMO_i;    
    //! LMO at next date.
    Data<T, 1> LMO_f;    
    //! LMO buffer.
    Data<T, 1> FileLMO_i;
    //! LMO buffer.
    Data<T, 1> FileLMO_f;

    //! WindDirectionInt at current date.
    Data<T, 1> WindDirectionInter_i;    
    //! WindDirectionInt at next date.
    Data<T, 1> WindDirectionInter_f;    
    //! WindDirectionInt buffer.
    Data<T, 1> FileWindDirectionInter_i;
    //! WindDirectionInt buffer.
    Data<T, 1> FileWindDirectionInter_f;

    //! WindSpeedInt at current date.
    Data<T, 1> WindSpeedInter_i;    
    //! WindSpeedInt at next date.
    Data<T, 1> WindSpeedInter_f;    
    //! WindSpeedInt buffer.
    Data<T, 1> FileWindSpeedInter_i;
    //! WindSpeedInt buffer.
    Data<T, 1> FileWindSpeedInter_f;
    
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

    //! PBLHInt at current date.
    Data<T, 1> PBLHInter_i;    
    //! PBLHInt at next date.
    Data<T, 1> PBLHInter_f;    
    //! PBLHInt buffer.
    Data<T, 1> FilePBLHInter_i;
    //! PBLHInt buffer.
    Data<T, 1> FilePBLHInter_f;

    //! USTInt at current date.
    Data<T, 1> USTInter_i;    
    //! USTInt at next date.
    Data<T, 1> USTInter_f;    
    //! USTInt buffer.
    Data<T, 1> FileUSTInter_i;
    //! USTInt buffer.
    Data<T, 1> FileUSTInter_f;

    //! LMOInt at current date.
    Data<T, 1> LMOInter_i;    
    //! LMOInt at next date.
    Data<T, 1> LMOInter_f;    
    //! LMOInt buffer.
    Data<T, 1> FileLMOInter_i;
    //! LMOInt buffer.
    Data<T, 1> FileLMOInter_f;
    
    /*** Species ***/

    //! Henry constants (mol / L / atm).
    map<string, T> henry_constant;
    //! Gas phase diffusivities (cm^2 / s).
    map<string, T> gas_phase_diffusivity;
    //! Constant scavenging coefficient (s^{-1}).
    map<string, T> scavenging_constant;

    //! Alpha scaling factor.
    map<string, T> alpha;
    //! Beta scaling factor..
    map<string, T> beta;
    
    //! Pointer to the current intersection.
    typename vector<Intersection<T>* >::iterator current_intersection;
    //! Intersections list.
    vector<Intersection<T>*> IntersectionVector;

    /*** Concentration ***/
    Data<T, 2> StreetConcentration;

    /*** New variables to particles ***/
    T min_water_drainage;
    string particles_dry_velocity_option;
    string brownian_diffusion_resistence_option;
    T max_rain;
    T f0_hdv; //s-1
    T f0_ldv; //s-1
    string file_traffic_flow;
    string file_traffic_speed;
    Array<T, 1> nb_hdv_hour;
    Array<T, 1> nb_ldv_hour;
    Array<T, 1> u_hdv;
    Array<T, 1> u_ldv;
    int Nh_week;
    //! diameter for each species.
    map<string, T> diameter;
    //! mass density of each species.
    map<string, T> mass_density;
    //! Drainage efficiency coefficient factor.
    map<string, T> drainage_efficiency;
    
    
    int rank;

  public:

    /*** Constructors and destructor ***/

    StreetNetworkTransport();
    StreetNetworkTransport(string config_file);
    virtual ~StreetNetworkTransport();

    /*** Initializations ***/

    void ReadConfiguration();
    void ReadStreetData();
    void CheckConfiguration();
    void DisplayConfiguration();
    
    void Init();
    void InitAllData();
    void InitStep();

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
    void SetStreetMeteoTransport(int street_index, T wind_direction, T wind_speed, 
				 T pblh, T ust, T lmo, T pressure, T temperature, T rain);
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
    void ComputeStreetConcentrationNoStationary();

    // Deposition
    void ComputeDryDepositionVelocities();
    string DepositionVelocityName(int s) const;
    int DepositionVelocityGlobalIndex(int s) const;
    void SetStreetDryDepositionFlux();
    void SetWallDryDepositionFlux();
    void SetStreetDryDepositionRate();
    void SetWallDryDepositionRate();
    // Scavenging
    void ComputeScavengingCoefficient();
    string ScavengingName(int s) const;
    int ScavengingGlobalIndex(int s) const;
    void SetStreetScavengingFlux();
    void SetStreetScavengingRate();
    bool HasScavenging(int s) const;
    bool HasScavenging(string name) const;
    bool ScavengingIndex(int s) const;
    bool ScavengingIndex(string name) const;
    
    void InitStep(T& sub_delta_t,
		  const T sub_delta_t_min,
		  const T transfer_velocity,
		  const T temp,
		  const T outgoing_flux,
		  const T street_volume,
		  const Array<T, 1> concentration_array,
		  const Array<T, 1> background_concentration_array,
		  const Array<T, 1> emission_rate_array,
		  const Array<T, 1> inflow_rate_array,
		  const Array<T, 1> deposition_flux_array,
		  const Array<T, 1> street_deposition_flux_array,
		  const Array<T, 1> scavenging_flux_array,
		  const T resuspension_factor,
		  const Array<T, 1> washoff_factor_array,
		  const Array<T, 1> street_deposited_mass_array,
		  const int Nsp,
		  const int bin);

    void ETRConcentration(const T transfer_velocity,
			  const T temp,
			  const T outgoing_flux,
			  const T street_volume,
			  const Array<T, 1> concentration_array,
			  Array<T, 1>& concentration_array_tmp, 
			  const Array<T, 1> background_concentration_array,
			  const Array<T, 1> emission_rate_array,
			  const Array<T, 1> inflow_rate_array,
			  const Array<T, 1> deposition_flux_array,
			  const Array<T, 1> street_deposition_flux_array,
			  const Array<T, 1> scavenging_flux_array,
			  const T resuspension_factor,
			  const Array<T, 1> washoff_factor_array,
			  const Array<T, 1> street_surface_deposited_mass_array,
			  Array<T, 1>& new_concentration_array,
			  const T sub_delta_t,
			  const int Nsp,
			  const int bin,
			  Array<T, 1>& new_street_surface_deposited_mass_array);

    void AdaptTimeStep(const Array<T, 1> new_concentration_array,
		       const Array<T, 1>concentration_array_tmp,
		       const T sub_delta_t_init,
		       const T sub_delta_t_min,
		       const T sub_delta_t_max,
		       T& sub_delta_t);

    void CalculDCDT(const T transfer_velocity,
		    const T temp,
		    const T outgoing_flux,
		    const T street_volume,
		    const Array<T, 1> concentration_array,
		    const Array<T, 1> background_concentration_array,
		    const Array<T, 1> emission_rate_array,
		    const Array<T, 1> inflow_rate_array,
		    const Array<T, 1> deposition_flux_array,
		    const Array<T, 1> street_deposition_flux_array,
		    const Array<T, 1> scavenging_flux_array,
		    const T resuspension_factor,
		    const Array<T, 1> washoff_factor_array,
		    const Array<T, 1> street_surface_deposited_mass,
		    Array<T, 1>& dcdt,
		    const int Nsp,
		    const int bin);

    void CalculDMDT(const Array<T, 1> street_surface_deposited_mass,
		    const Array<T, 1> concentration_array,
		    const Array<T, 1> street_deposition_flux_array,
		    const Array<T, 1> street_scavenging_flux_array,
		    const T resuspension_factor,
		    const Array<T, 1> washoff_factor_array,
		    Array<T, 1>& dmdt,
		    const int Nsp,
		    const int bin);

    void RosenbrockConcentration(const T transfer_velocity,
				 const T temp,
				 const T outgoing_flux,
				 const T street_volume,
				 const Array<T, 1> concentration_array,
				 Array<T, 1>& concentration_array_tmp,
				 const Array<T, 1> background_concentration_array,
				 const Array<T, 1> emission_rate_array,
				 const Array<T, 1> inflow_rate_array,
				 const Array<T, 1> deposition_flux_array,
				 Array<T, 1>& new_concentration_array,
				 const T sub_delta_t,
				 const T inflow_flux,
				 const Array<T, 1> street_deposition_flux_array,
				 const Array<T, 1> scavenging_flux_array,
				 const T resuspension_factor,
				 const Array<T, 1> washoff_factor_array,
				 const Array<T, 1> initial_deposited_mass,
				 Array<T, 1>& new_street_surface_deposited_mass_array,
				 const int Nsp,
				 const int bin);
    
    void Jacobian(const T inflow_flux,
		  const T outgoing_flux,
		  const T temp,
		  const Array<T, 1> deposition_flux_array,
		  Array<T, 1>& J,
		  const T street_volume,
		  int Nsp,
		  int bin);
    
    void Jacobian(const T resuspension_factor,
		  const Array<T, 1> washoff_factor_array,
		  const Array<T, 1> deposition_flux_array,
		  Array<T, 1>& J,
		  const T street_volume,
		  int Nsp,
		  int bin);

    void CalculWashoffFactor();
    void CalculAerosolResuspensionFactor();
    void ComputeStreetSurfaceDepositedMass();
    void IsStationary(bool& is_stationary);
    void OutputSaver();
    void InitOutputSaver();
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKTRANSPORTINERTPARTICLE_HXX
#endif
