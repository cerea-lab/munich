#ifndef POLYPHEMUS_FILE_MODELS_STREETNETWORKAEROSOL_HXX

#include "StreetNetworkChemistry.cxx"
#include "StreetAerosol.cxx"

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

  //////////////////
  // STREET //
  //////////////////


  //! StreetNetworkAerosol model.
  // template<class T, class ClassChemistry>
  // class StreetNetworkAerosol: public StreetNetworkChemistry<T, ClassChemistry>,
  //                             public BaseModuleParallel
  template<class T, class ClassChemistry>
  class StreetNetworkAerosol: public StreetNetworkChemistry<T, ClassChemistry>
  {

  protected:

    static const T pi;

    //Number//
    string number_computation_option;
    
    //***** Chemistry *****//
    ClassChemistry Chemistry_;
    string option_chemistry;
    T Rho_tmp_aer;
    T rho_mean_aer;
    vector<T> Rho_species_aer; // g.cm-3 -> 10-6 ug.um-3 
    
    //***** Size sections discretization *****//
    vector<string> bin_list;
    Array<T, 1> BinBound_aer;
    
    //***** Mass fraction discretization *****//
    vector<string> fraction_list;
    int Nfraction_aer;
    Array<T, 1> Fractionbound_aer;
    T fixed_density_aer;    

    //***** Initial conditions *****//
    
    //!Tag of input data format for initial conditions
    string ic_format;
    //! Number of aerosol species with initial conditions.
    int Ns_ic_aer;
    //! Number of aerosol number bins with initial conditions.
    int Nb_ic_aer;    
    //! List of aerosol species (with their bins) with initial conditions.
    vector<pair<string, vector<int> > > species_list_ic_aer;
    //! List of aerosol number bins with initial conditions.
    vector<int> ic_bin_list_aer;
    
    /*** Deposision ***/
    T max_rain;
    string washoff_option;
    string particles_dry_velocity_option;
    string brownian_diffusion_resistence_option;
    T min_water_drainage;
    //vector<int> bin_list_dep_aer;
    int Nb_dep_aer; //number of size sections deposed
    int Nbin_dep_aer; //number of size sections deposed x number of particles compositions
    //! Dry deposition fluxes at current date.
    Data<T, 2> StreetDryDepositionFlux_aer;
    //! Dry wall deposition fluxes at current date.
    Data<T, 2> WallDryDepositionFlux_aer;
    //! Dry tree leaves deposition fluxes at current date.
    Data<T, 2> TreeDryDepositionFlux_aer;
    //! Dry deposition rate at current date.
    Data<T, 2> StreetDryDepositionRate_aer;
    //! Dry wall deposition rate at current date.
    Data<T, 2> WallDryDepositionRate_aer;
    //! Dry tree leaves deposition rate at current date.
    Data<T, 2> TreeDryDepositionRate_aer;
    //! List of species with deposition.
    vector<pair<string, vector<int> > > species_list_dep_aer;
    //! List of aerosol bins with deposition.
    vector<int> dep_bin_list_aer;
    
    //! Dry wall  deposition rate at current date.
    Data<T, 3> StreetSurfaceDepositedMass_aer;
    //! Dry wall  deposition rate at current date.
    Data<T, 3> StreetResuspensionRate_aer;
    //! Dry wall  deposition rate at current date.
    Data<T, 3> StreetSurfaceDryDepositionRate_aer;
    //! Dry wall  deposition rate at current date.
    Data<T, 3> StreetWashoffRate_aer;
    

    //! Dry wall deposition rate at current date.
    //Data<T, 3> InitialStreetSurfaceDepositedMass_aer;

    //! Drainage efficiency coefficient factor.
    map<string, T> drainage_efficiency;
    
    /*** Scavenging ***/
    //! Dissolution heat (kcal / mol at 298K).
    map<string, T> dissolution_heat;
    //! Number of aerosol bins with scavenging.
    int Nb_scav_aer;
    int Nbin_scav_aer;
    //! Scavenging coefficient at current date.
    Data<T, 2> ScavengingCoefficient_aer;
    //! Scavenging fluxes at current date.
    Data<T, 2> StreetScavengingFlux_aer;
    //! Scavenging rate at current date.
    Data<T, 2> StreetScavengingRate_aer;
    //! List of aerosol bins with scavenging.
    vector<int> scav_bin_list_aer;

    /*** Resuspension and Non-exhaust emissions***/

    string file_resuspension;
    T f0_hdv; //veh-1
    T f0_lcv; //veh-1
    T f0_2R; //veh-1
    T f0_pc; //veh-1
    T mean_speed_2R; //veh.h-1
    T mean_speed_HDV; //veh.h-1
    T mean_speed_PC; //veh.h-1
    T mean_speed_LCV; //veh.h-1
    T mean_speed_highway_2R; //veh.h-1
    T mean_speed_highway_HDV; //veh.h-1
    T mean_speed_highway_PC; //veh.h-1
    T mean_speed_highway_LCV; //veh.h-1

    string wet_diameter_option;
    
    /*** Grids ***/
    RegularGrid<T> GridS3D_aer;
    RegularGrid<T> GridB3D_aer;
    
    RegularGrid<T> GridS2D_aer;
    RegularGrid<T> GridB2D_aer;
    
    RegularGrid<T> GridS_emis_aer;
    RegularGrid<T> GridB_emis_aer; 
    
    RegularGrid<T> GridS_bg_aer;
    RegularGrid<T> GridB_bg_aer;     

    RegularGrid<T> GridS_dep_aer;
    RegularGrid<T> GridB_dep_aer;   

    RegularGrid<T> GridS_scav_aer;
    RegularGrid<T> GridB_scav_aer;

    RegularGrid<T> GridB2D_aer_i;
    RegularGrid<T> GridB3D_aer_i;
    
    //! 2D grid for groups.
    RegularGrid<T> GridG2D_aer;

    //! 3D grid for streets.
    RegularGrid<T> GridST3D;
    
    
    /*** Concentration ***/
    Data<T, 3> StreetConcentration_aer;  
    Data<T, 2> StreetNumberConcentration;

    Data<T, 3> StreetConcentration_aer_i;  
    Data<T, 2> StreetNumberConcentration_i; 
    
    /*** Background ***/
    //! Aerosol background concentration at current date.
    Data<T, 3> Background_aer_i;    
    //! Aerosol emissions at next date.
    Data<T, 3> Background_aer_f;    
    //! Aerosol emissions buffer.
    Data<T, 3> FileBackground_aer_i;
    //! Aerosol emissions buffer.
    Data<T, 3> FileBackground_aer_f;
    //! Number emissions at current date for aerosols.
    Data<T, 2> NumberBackground_aer_i;
    //! Number emissions at next date for aerosols.
    Data<T, 2> NumberBackground_aer_f;
    //! Number emissions buffer.
    Data<T, 2> FileNumberBackground_aer_i;
    //! Number emissions buffer.
    Data<T, 2> FileNumberBackground_aer_f;      
    
    
    int Nt_traffic;
    string aerosol_bg_format;
    int Nt_bg_aer;
    vector<pair<string, vector<int> > > species_list_bg_aer;
    vector<int> bg_bin_list_aer;
    int Ns_bg_aer;
    int Nb_bg_aer; 
    Array<T, 4> bg_concentration_aer;
    Array<T, 3> total_bin_bg_aer;    
    Array<T, 3> bg_number_concentration;   
    
    /*** Emission ***/
    //! Aerosol emissions at current date.
    Data<T, 3> Emission_aer_i;    
    //! Aerosol emissions at next date.
    Data<T, 3> Emission_aer_f;    
    //! Aerosol emissions buffer.
    Data<T, 3> FileEmission_aer_i;
    //! Aerosol emissions buffer.
    Data<T, 3> FileEmission_aer_f;
    //! Number emissions at current date for aerosols.
    Data<T, 2> NumberEmission_aer_i;
    //! Number emissions at next date for aerosols.
    Data<T, 2> NumberEmission_aer_f;
    //! Number emissions buffer.
    Data<T, 2> FileNumberEmission_aer_i;
    //! Number emissions buffer.
    Data<T, 2> FileNumberEmission_aer_f;      
    
    string aerosol_emission_format;
    int Nt_emis_aer;
    vector<pair<string, vector<int> > > species_list_emis_aer;
    int Ns_emis_aer;
    vector<int> emis_bin_list_aer;
    int Nb_emis_aer;
    Array<T, 4> emission_aer;
    Array<T, 3> total_bin_emission_aer;
    Array<T, 3> number_emission_aer;  
    
    /*** New meteo data ***/
    //! LiquidWaterContent at current date.
    Data<T, 1> LiquidWaterContent_i;    
    //! LiquidWaterContent at next date.
    Data<T, 1> LiquidWaterContent_f;    
    //! LiquidWaterContent buffer.
    Data<T, 1> FileLiquidWaterContent_i;
    //! LiquidWaterContent buffer.
    Data<T, 1> FileLiquidWaterContent_f;
    //! RelativeHumidity at current date.    
    Data<T, 1> RelativeHumidity_i;    
    //! RelativeHumidity at next date.
    Data<T, 1> RelativeHumidity_f;

    /*** Traffic data ***/
    //! 2R road traffic at current date.
    Data<T, 1> RoadTraffic_2R_i;    
    //! 2R road traffic at next date.
    Data<T, 1> RoadTraffic_2R_f;
    //! 2R road traffic buffer.
    Data<T, 1> FileRoadTraffic_2R_i;
    //! 2R road traffic buffer.
    Data<T, 1> FileRoadTraffic_2R_f;

    //! HDV road traffic at current date.
    Data<T, 1> RoadTraffic_HDV_i;    
    //! HDV road traffic at next date.
    Data<T, 1> RoadTraffic_HDV_f;
    //! HDV road traffic buffer.
    Data<T, 1> FileRoadTraffic_HDV_i;
    //! HDV road traffic buffer.
    Data<T, 1> FileRoadTraffic_HDV_f;

    //! PC road traffic at current date.
    Data<T, 1> RoadTraffic_PC_i;    
    //! PC road traffic at next date.
    Data<T, 1> RoadTraffic_PC_f;
    //! PC road traffic buffer.
    Data<T, 1> FileRoadTraffic_PC_i;
    //! PC road traffic buffer.
    Data<T, 1> FileRoadTraffic_PC_f;

    //! LCV road traffic at current date.
    Data<T, 1> RoadTraffic_LCV_i;    
    //! LCV road traffic at next date.
    Data<T, 1> RoadTraffic_LCV_f;
    //! LCV road traffic buffer.
    Data<T, 1> FileRoadTraffic_LCV_i;
    //! LCV road traffic buffer.
    Data<T, 1> FileRoadTraffic_LCV_f;
    
    //OLD 
    // Array<T, 2> rain;
    // Array<T, 2> relative_humidity;  
    // Array<T, 2> liquidwatercontent;
    Data<T, 1> pH;

    // //! Mass density (\mu g.\mu m^-3).
    // vector<T> Mass_Density_aer;
    
    //!#### SSH #### 

    //! Number of non-organic aerosol species.
    int n_nonorganic;

    //! Number of organic aerosol species.
    int n_organics;    
    
    int ssh_Ns_aer_nolayer;
    
    //! Number of aerosol layers
    int N_layer;
    //! i_hydrophilic option from ssh
    int i_hydrophilic;

    Array<string, 1> species_list_aer_nolayer;    

    Array<T, 1> ssh_mass_density_layers;
    
    Array<int, 1> index_species_layers;

    Array<int, 1> aerosol_species_group_relation;

    Array<T, 3> composition_bounds;

    Array<T, 1> aerosol_type;
    
    
  public:

    /*** Constructors and destructor ***/

    StreetNetworkAerosol();
    StreetNetworkAerosol(string config_file);
    virtual ~StreetNetworkAerosol();  

    /*** Configuration ***/
    //! Threshold on liquid water content for clouds (g / m^3).
    T lwc_cloud_threshold;

    /*** Initializations ***/

    void ReadConfiguration();
    void CheckConfiguration();
    void DisplayConfiguration();
    void Allocate();    

    void Init();
    void InitStreet();
    void InitStep();
    void InitAllData();
//     void InitData(string input_file, Array<T, 2>& input_data);

    void OutputSaver();  
    
    /*** Access Methods ***/

    // void SetStreetAdditionalMeteoAerosol(int street_index, T attenuation,
    // 					 T specific_humidity, T pressure,
    // 					 T temperature, T rain);
    
    Data<T, 3>& GetStreetConcentration_aer();
    Data<T, 2>& GetStreetNumberConcentration();
    void ComputeStreetSurfaceDepositedMass_aer();
    void SetStreetOutput_aer();
    void SetStreetDryDepositionFlux_aer();
    void SetWallDryDepositionFlux_aer();
    void SetTreeDryDepositionFlux_aer();
    void SetStreetDryDepositionRate_aer();
    void SetWallDryDepositionRate_aer();
    void SetTreeDryDepositionRate_aer();
    void SetStreetScavengingFlux_aer();
    void SetStreetScavengingRate_aer();
    void SetTotalDryDepositionRate_aer();
    void SetStreetScavengingFluxOverCanopy_aer(int street_index, T scavenging_flux_overcanopy, int b);
    T GetStreetScavengingFlux_aer(int street_index, int b);
    Data<T, 2>& GetStreetScavengingFlux_aer();
    Data<T, 2>& GetStreetScavengingRate_aer();
    void SetpH();

    void SetStreetSurfaceDepositedMass_aer();
    void SetStreetResuspensionRate_aer();
    void SetStreetWashoffRate_aer();
    void SetStreetSurfaceDryDepositionRate_aer();
    
    
    /*** Computational Methods ***/
    
    void Forward();
    void InitMeteo();
    void ComputeMassBalance();
    void ComputeInflowRateExtended();
    void ComputeStreetConcentration();
    void ComputeStreetConcentrationNoStationary();
    void AdaptTimeStep(const Array<T, 1> new_concentration_array,
		       const Array<T, 1>concentration_array_tmp,
		       const Array<T, 2> new_concentration_array_aer,
		       const Array<T, 2>concentration_array_tmp_aer,
		       const T sub_delta_t_init,
		       const T sub_delta_t_min,
		       const T sub_delta_t_max,
		       T& sub_delta_t);
    
    void AdaptTimeStep(const Array<T, 1> new_concentration_array,
		       const Array<T, 1>concentration_array_tmp,
		       const Array<T, 2> new_concentration_array_aer,
		       const Array<T, 2>concentration_array_tmp_aer,
		       const Array<T, 1> new_number_concentration_array,
		       const Array<T, 1>number_concentration_array_tmp,
		       const T sub_delta_t_init,
		       const T sub_delta_t_min,
		       const T sub_delta_t_max,
		       T& sub_delta_t);
    void InitWetDiameter_aer(Array<T, 1>& WetDiameter_aer_);    
    void ComputeDryDepositionVelocities_aer();
    //void ComputeScavengingCoefficient();
    void ComputeScavengingCoefficient_aer();

    /*** Access Methods: Chemistry ***/
    void Chemistry();

    void Chemistry(Date current_date_tmp,
		   T sub_delta_t,
		   T attenuation_,
		   T specific_humidity_,
		   T pressure_,
		   T temperature_,
		   T longitude_,
		   T latitude_,
		   T rain_,
		   T liquidwatercontent_,
		   Array<T, 1> photolysis_rate,
		   Array<T, 1>& concentration_array,
		   Array<T, 2>& concentration_array_aer,
		   Array<T, 1>& number_concentration,
		   Array<T, 1> wet_diameter_aer,
		   Array<T, 1>& wet_diameter_aer_loc);

    T ComputeDensity(Data <T, 1> Conc_aer_tmp,vector<T> Rho_species, T TotalMass, int Ns);
    void CalculNumberResuspension();
    void ComputeNumberEmission_aer(int b);
    void ComputeNumberBackground_aer(int b);
    int Bin_to_size_index_aer(int b) const; 
    int NumberEmissionIndex_aer(int x) const;   
    bool HasNumberEmission_aer(int x) const;
    int Bin_index_translate_aer(int s, int b) const;
    int Bin_to_composition_index(int b) const;
    int NumberBackgroundIndex_aer(int x) const;
    bool HasNumberBackground_aer(int x) const;    

    void SetStreetConcentration_aer();
    void SetStreetNumberConcentration_aer();
    void InitChemistry();
    void InitInflowRate();
    void SetInitialStreetConcentration();
    int FindExternalCompositionID(string species);
    
    bool HasNumberConcentration_aer();
    bool HasScavenging_aer(int x) const;
    int ScavengingIndex_aer(int b) const;

    int DepositionVelocityGlobalIndex_aer(int s) const;
    string DepositionVelocityName_aer(int s) const;
    void CalculWashoffFactor();
    /* Resuspension */
    void CalculAerosolResuspensionFactor();
    void SetStreetRoadTraffic(Data<T, 1> RoadTraffic_2R_f,
			      Data<T, 1> RoadTraffic_HDV_f,
			      Data<T, 1> RoadTraffic_PC_f,
			      Data<T, 1> RoadTraffic_LCV_f);
    void AppendBinList_aer(int b, vector<int>& bin_list_aer); // YK
    void InitStreetConc();
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKAEROSOL_HXX
#endif
