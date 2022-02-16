#ifndef POLYPHEMUS_FILE_MODELS_STREETNETWORKTRANSPORT_CXX

//////////////
// INCLUDES //

#include "StreetNetworkTransport.hxx"

// INCLUDES //
//////////////
  
namespace Polyphemus
{

  template<class T>
  const T StreetNetworkTransport<T>::pi = acos(-1);

  template<class T>
  const T StreetNetworkTransport<T>::earth_radius = 6371229.;

  template<class T>
  const T StreetNetworkTransport<T>::karman = 0.41;

  template<class T>
  const T StreetNetworkTransport<T>::nu = 0.15;

  template<class T>
  const T StreetNetworkTransport<T>::Pr = 0.74;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds the model. Nothing else is performed.
  */
  template<class T>
  StreetNetworkTransport<T>::StreetNetworkTransport():
    BaseModel<T>()
  {
  }


  //! Main constructor.
  /*!
    \param config_file configuration filename.
  */
  template<class T>
  StreetNetworkTransport<T>::StreetNetworkTransport(string config_file):
    BaseModel<T>(config_file)
  {
    this->D2_map["StreetConcentration"] = &StreetConcentration;
    this->D2_map["StreetDryDepositionVelocity"] =
      &StreetDryDepositionVelocity;
  }
  

  //! Destructor.
  template<class T>
  StreetNetworkTransport<T>::~StreetNetworkTransport()
  {
    ClearStreetVector();
    ClearIntersectionVector();
  }

  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T>
  void StreetNetworkTransport<T>::ReadConfiguration()
  {

    /*** Display ***/
    this->config.SetSection("[display]");
    // Should the configuration be displayed on screen?
    if (this->config.Check("Show_configuration"))
      this->config.PeekValue("Show_configuration",
                             this->option_process["show_configuration"]);
    else
      this->option_process["show_configuration"] = true;
      
    /*** Dates ***/
    this->config.SetSection("[domain]");
    this->Date_min = this->config.PeekValue("Date_min");
    this->config.PeekValue("Delta_t", "> 0", this->Delta_t);
    this->config.PeekValue("Nt", "positive", this->Nt);

    // /*** Species ***/
    this->config.PeekValue("Species", this->file_species);
    //! Opens the file that describes species.
    ConfigStream species_stream(this->file_species);
    //! Section "[species]" contains all species names.
    species_stream.SetSection("[species]");
    while (!species_stream.IsEmpty())
      this->species_list.push_back(species_stream.GetElement());
    this->Ns = int(this->species_list.size());

    /*** Options ***/
    
    this->config.SetSection("[options]");
    if (this->config.Check("With_transport"))
      this->config.PeekValue("With_transport",
                             this->option_process["with_transport"]);
    else
      this->option_process["with_transport"] = true;
    if (this->config.Check("With_deposition"))
      this->config.PeekValue("With_deposition",
                             this->option_process["with_deposition"]);
    else
      this->option_process["with_deposition"] = true;

    if (this->option_process["with_deposition"])
      {
        if (this->config.Check("Collect_dry_flux"))
          this->config.PeekValue("Collect_dry_flux",
                                 this->option_process["collect_total_dry_flux"]);
        else
          this->option_process["collect_total_dry_flux"] = false;
      }

    if (this->config.Check("With_scavenging"))
      this->config.PeekValue("With_scavenging",
                             this->option_process["with_scavenging"]);
    else
      this->option_process["with_scavenging"] = true;
    if (this->option_process["with_scavenging"])
      this->scavenging_model = "microphysical";
    else
      this->scavenging_model = "none";
    this->scavenging_model = lower_case(this->scavenging_model);

    if (this->config.Check("With_initial_condition"))
      this->config.PeekValue("With_initial_condition",
                             this->option_process["with_initial_condition"]);
    else
      this->option_process["with_initial_condition"] = false;

    this->config.SetSection("[street]");
    
    this->config.PeekValue("Transfer_parameterization",
                           "Sirane | Schulte | Wang", option_transfer);
    this->config.PeekValue("Mean_wind_speed_parameterization",
                           "Sirane | Exponential | Wang", option_ustreet);
    this->config.PeekValue("Building_height_wind_speed_parameterization",
                           "Sirane | Macdonald", option_uH);
    this->config.PeekValue("Compute_Macdonald_from",
                           "Uref | Ustar", option_uH_macdonald);

    if (this->config.Check("Zref"))
      this->config.PeekValue("Zref", "> 0", zref);
    else
      zref = 30.0;

      if (this->config.Check("z0_surface"))
      this->config.PeekValue("z0_surface", "> 0", z0s);
    else
      z0s = 0.05;
   
    this->config.PeekValue("With_horizontal_fluctuation",
			   this->option_process["with_horizontal_fluctuation"]);
    if (this->config.Check("Minimum_Street_Wind_Speed"))
      this->config.PeekValue("Minimum_Street_Wind_Speed",
                             ">= 0.1", ustreet_min);
    else
      ustreet_min = 0.1;

    if (this->config.Check("With_local_data"))
      this->config.PeekValue("With_local_data",
                             this->option_process["with_local_data"]);
    else
      this->option_process["with_local_data"] = true;

    if (this->option_process["with_deposition"] or
        option_uH == "Macdonald")
      {
        if (this->config.Check("Building_density")) 
          this->config.PeekValue("Building_density",
                                 this->building_density);
        else
          this->building_density = 0.4;
        
      }
    
    this->config.PeekValue("With_stationary_hypothesis",
			   this->option_process["with_stationary_hypothesis"]);

    if(!this->option_process["with_stationary_hypothesis"])
      {
	
        this->config.PeekValue("Numerical_method_parameterization",
                               "ETR | Rosenbrock", option_method);
        if (this->config.Check("sub_delta_t_min"))
          this->config.PeekValue("sub_delta_t_min", ">= 0.1",
                                 sub_delta_t_min);
        else
          sub_delta_t_min = 1.0;
      }

    this->config.PeekValue("With_tree_aerodynamic",
                           this->option_process["with_tree_aerodynamic"]);
    //! Aerodynamic effect of trees must be computed with Wang parameterization
    if (this->option_process["with_tree_aerodynamic"] and (option_transfer != "Wang" or option_ustreet != "Wang"))
      throw ("Choose option Wang for Mean_wind_speed_parameterization and Transfer_parameterization if option with_tree_aerodynamic is activated");

    this->config.PeekValue("With_tree_deposition",
                           this->option_process["with_tree_deposition"]);
    if (this->option_process["with_tree_deposition"] and !this->option_process["with_deposition"])
      throw ("Option With_deposition has to be activated to compute dry deposition on tree leaves.");

    if (this->option_process["with_tree_deposition"])
      this->config.PeekValue("Cuticular_resistance_option",
                             "zhang | wesely", option_Rcut);

    if (this->option_process["with_tree_aerodynamic"] or this->option_process["with_tree_deposition"])
      {
        if (this->config.Check("Tree_drag_coefficient"))
          this->config.PeekValue("Tree_drag_coefficient", "> 0", Cdt);
        else
          Cdt = 0.2;
      }

    // Tree type (LUC) to compute deposition on tree leaves
    if (this->option_process["with_tree_deposition"])
      {
        this->config.PeekValue("Tree_LUC", tree_LUC);
        file_tree_deposition = "tree_deposition_param.dat";
      }

    /*** Input files ***/

    //! The configuration-file path is the field "Data_description" in the main
    //! configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    //! Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    //! Meteorological files.
    if (this->option_process["with_local_data"])
      {
        this->input_files["meteo"].Read(data_description_file, "meteo");
        data_description_stream.SetSection("[meteo]");
        data_description_stream.PeekValue("Nt", "> 0", Nt_meteo);
      }

    //! Emission files.
    this->input_files["emission"].Read(data_description_file, "emission");
    data_description_stream.SetSection("[emission]");
    data_description_stream.PeekValue("Nt", "> 0", Nt_emis);
    for (map<string, string>::iterator i
           = this->input_files["emission"].Begin();
         i != this->input_files["emission"].End(); i++)
      species_list_emis.push_back(i->first);
    Ns_emis = int(species_list_emis.size());

    //! Background concentrations files.
    if (this->option_process["with_local_data"])
      {
        this->input_files["background_concentration"].Read(data_description_file,
                                                           "background_concentration");
        data_description_stream.SetSection("[background_concentration]");
        data_description_stream.PeekValue("Nt", "> 0", Nt_background);
        for (map<string, string>::iterator i
               = this->input_files["background_concentration"].Begin();
             i != this->input_files["background_concentration"].End(); i++)
          species_list_background.push_back(i->first);
        Ns_background = int(species_list_background.size());
      }

    //! Initial conditions
    if (this->option_process["with_initial_condition"])
      this->input_files["initial_condition"].ReadFiles(data_description_file,
                                                       "initial_condition");
    else
      this->input_files["initial_condition"].Empty();      
    //        data_description_stream.SetSection("[initial_conditions]");
    for (map<string, string>::iterator i
           = this->input_files["initial_condition"].Begin();
         i != this->input_files["initial_condition"].End(); i++)
      species_list_ic.push_back(i->first);
    Ns_ic = int(species_list_ic.size());

    // Deposition
    if (this->option_process["with_deposition"])
      this->input_files["deposition"]
        .ReadFields(data_description_file, "deposition");
    else
      this->input_files["deposition"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["deposition"].Begin();
	 i != this->input_files["deposition"].End(); i++)
      species_list_dep.push_back(i->first);
    Ns_dep = int(species_list_dep.size());

   // Scavenging.
    this->input_files["scavenging_coefficient"].Empty(); 
    if (this->option_process["with_scavenging"])
      {
        if (this->scavenging_model != "none")
          this->input_files["scavenging_coefficient"]
            .ReadFields(data_description_file, "scavenging");
    
        for (map<string, string>::iterator i
               = this->input_files["scavenging_coefficient"].Begin();
             i != this->input_files["scavenging_coefficient"].End(); i++)
          species_list_scav.push_back(i->first);
        Ns_scav = int(species_list_scav.size());
      }
    else
      Ns_scav = 0;
    
    /*** Species data ***/

    ConfigStream species_data_stream(this->GetSpeciesFile());
    string species;

    //Molecular weight []
    species_data_stream.SetSection("[molecular_weight]");
    while (!species_data_stream.IsEmpty())
      {
        species = species_data_stream.GetElement();
        species_data_stream.GetNumber(molecular_weight[species]);
      }

    //Reactivity []
    species_data_stream.SetSection("[reactivity]");
    while (!species_data_stream.IsEmpty())
      {
        species = species_data_stream.GetElement();
        species_data_stream.GetNumber(reactivity[species]);
      }


    // Alpha scaling factor for deposition.
    species_data_stream.SetSection("[alpha]");
    while (!species_data_stream.IsEmpty())
      {
        species = species_data_stream.GetElement();
        species_data_stream.GetNumber(alpha[species]);
      }

    // Beta scaling factor for deposition.
    species_data_stream.SetSection("[beta]");
    while (!species_data_stream.IsEmpty())
      {
        species = species_data_stream.GetElement();
        species_data_stream.GetNumber(beta[species]);
      }
	
    // Henry constants in mol / L / atm.
    species_data_stream.SetSection("[henry]");
    while (!species_data_stream.IsEmpty())
      {
        species = species_data_stream.GetElement();
        species_data_stream.GetNumber(henry_constant[species]);
      }
      
    // Gas phase diffusivity constants in cm^2 / s.
    species_data_stream.SetSection("[diffusivity]");
    while (!species_data_stream.IsEmpty())
      {
        species = species_data_stream.GetElement();
        species_data_stream.GetNumber(gas_phase_diffusivity[species]);
      }

    ReadStreetData();

    if (this->option_process["with_tree_aerodynamic"] or this->option_process["with_tree_deposition"])
      ReadTreeData();

    CheckConfiguration();

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif
    if (this->option_process["show_configuration"] && rank == 0)
      DisplayConfiguration();   
  }

  //! Display the configuration
  template<class T>
  void StreetNetworkTransport<T>::DisplayConfiguration()
  {

    cout << "Number of streets: " << total_nstreet << endl;
    cout << "Transfer_parameterization: " << option_transfer << endl;
    cout << "Mean_wind_speed_parameterization: " << option_ustreet << endl;
    cout << "Scavenging model: " << this->scavenging_model << endl;
    
    if(!this->option_process["with_stationary_hypothesis"])
      {
	cout << "Non-stationary approach to calculate concentrations in the streets." << endl;
	cout << "Numerical_method_parameterization: " << option_method << endl;
	cout << "sub_delta_t_min: " << sub_delta_t_min << endl;
      }
    else
      cout << "Stationary approach to calculate concentrations in the streets." << endl;

    if(this->option_process["with_deposition"])
      cout<<"With dry deposition of gas-phase pollutants in streets."<<endl;
    if(this->option_process["with_scavenging"])
      cout<<"With wet deposition of gas-phase pollutants in streets."<<endl;

    if(this->option_process["with_tree_aerodynamic"])
      cout<<"With tree aerodynamic effect in streets."<<endl;
    if(this->option_process["with_tree_deposition"])
      cout<<"With dry deposition of gas-phase pollutants and aerosols on tree leaves."<<endl;
  }
  
  //! Allocates memory.
  /*! Allocates the grids and the concentration Data.
   */
  template<class T>
  void StreetNetworkTransport<T>::Allocate()
  {
    
    GridST1D = RegularGrid<T>(total_nstreet);
    GridST2D = RegularGrid<T>(total_nstreet);
    
    GridINT2D = RegularGrid<T>(nintersection);
    GridS2D = RegularGrid<T>(this->Ns);

    /*** Meteo ***/
    if (this->option_process["with_local_data"])
      {
	Temperature_i.Resize(GridST2D);
	Temperature_f.Resize(GridST2D);
	FileTemperature_i.Resize(GridST2D);
	FileTemperature_f.Resize(GridST2D);
	
	Temperature_i.SetZero();
	Temperature_f.SetZero();
	FileTemperature_i.SetZero();
	FileTemperature_f.SetZero();

	Pressure_i.Resize(GridST2D);
	Pressure_f.Resize(GridST2D);
	FilePressure_i.Resize(GridST2D);
	FilePressure_f.Resize(GridST2D);
	
	Pressure_i.SetZero();
	Pressure_f.SetZero();
	FilePressure_i.SetZero();
	FilePressure_f.SetZero();

	Rain_i.Resize(GridST2D);
	Rain_f.Resize(GridST2D);
	FileRain_i.Resize(GridST2D);
	FileRain_f.Resize(GridST2D);
	
	Rain_i.SetZero();
	Rain_f.SetZero();
	FileRain_i.SetZero();
	FileRain_f.SetZero();

	WindDirection_i.Resize(GridST2D);
	WindDirection_f.Resize(GridST2D);
	FileWindDirection_i.Resize(GridST2D);
	FileWindDirection_f.Resize(GridST2D);
	
	WindDirection_i.SetZero();
	WindDirection_f.SetZero();
	FileWindDirection_i.SetZero();
	FileWindDirection_f.SetZero();

	WindSpeed_i.Resize(GridST2D);
	WindSpeed_f.Resize(GridST2D);
	FileWindSpeed_i.Resize(GridST2D);
	FileWindSpeed_f.Resize(GridST2D);
	
	WindSpeed_i.SetZero();
	WindSpeed_f.SetZero();
	FileWindSpeed_i.SetZero();
	FileWindSpeed_f.SetZero();

	PBLH_i.Resize(GridST2D);
	PBLH_f.Resize(GridST2D);
	FilePBLH_i.Resize(GridST2D);
	FilePBLH_f.Resize(GridST2D);
	
	PBLH_i.SetZero();
	PBLH_f.SetZero();
	FilePBLH_i.SetZero();
	FilePBLH_f.SetZero();

	UST_i.Resize(GridST2D);
	UST_f.Resize(GridST2D);
	FileUST_i.Resize(GridST2D);
	FileUST_f.Resize(GridST2D);
	
	UST_i.SetZero();
	UST_f.SetZero();
	FileUST_i.SetZero();
	FileUST_f.SetZero();

	LMO_i.Resize(GridST2D);
	LMO_f.Resize(GridST2D);
	FileLMO_i.Resize(GridST2D);
	FileLMO_f.Resize(GridST2D);
	
	LMO_i.SetZero();
	LMO_f.SetZero();
	FileLMO_i.SetZero();
	FileLMO_f.SetZero();

	WindDirectionInter_i.Resize(GridINT2D);
	WindDirectionInter_f.Resize(GridINT2D);
	FileWindDirectionInter_i.Resize(GridINT2D);
	FileWindDirectionInter_f.Resize(GridINT2D);
	
	WindDirectionInter_i.SetZero();
	WindDirectionInter_f.SetZero();
	FileWindDirectionInter_i.SetZero();
	FileWindDirectionInter_f.SetZero();

	WindSpeedInter_i.Resize(GridINT2D);
	WindSpeedInter_f.Resize(GridINT2D);
	FileWindSpeedInter_i.Resize(GridINT2D);
	FileWindSpeedInter_f.Resize(GridINT2D);
	
	WindSpeedInter_i.SetZero();
	WindSpeedInter_f.SetZero();
	FileWindSpeedInter_i.SetZero();
	FileWindSpeedInter_f.SetZero();

	PBLHInter_i.Resize(GridINT2D);
	PBLHInter_f.Resize(GridINT2D);
	FilePBLHInter_i.Resize(GridINT2D);
	FilePBLHInter_f.Resize(GridINT2D);
	
	PBLHInter_i.SetZero();
	PBLHInter_f.SetZero();
	FilePBLHInter_i.SetZero();
	FilePBLHInter_f.SetZero();

	USTInter_i.Resize(GridINT2D);
	USTInter_f.Resize(GridINT2D);
	FileUSTInter_i.Resize(GridINT2D);
	FileUSTInter_f.Resize(GridINT2D);
	
	USTInter_i.SetZero();
	USTInter_f.SetZero();
	FileUSTInter_i.SetZero();
	FileUSTInter_f.SetZero();

	LMOInter_i.Resize(GridINT2D);
	LMOInter_f.Resize(GridINT2D);
	FileLMOInter_i.Resize(GridINT2D);
	FileLMOInter_f.Resize(GridINT2D);
	
	LMOInter_i.SetZero();
	LMOInter_f.SetZero();
	FileLMOInter_i.SetZero();
	FileLMOInter_f.SetZero();

	SpecificHumidity_i.Resize(GridST2D);
	SpecificHumidity_f.Resize(GridST2D);
	FileSpecificHumidity_i.Resize(GridST2D);
	FileSpecificHumidity_f.Resize(GridST2D);
	
	SpecificHumidity_i.SetZero();
	SpecificHumidity_f.SetZero();
	FileSpecificHumidity_i.SetZero();
	FileSpecificHumidity_f.SetZero();

        // New meteo for gas deposition on tree leaves
        SolarRadiation_i.Resize(GridST2D);
        SolarRadiation_f.Resize(GridST2D);
        FileSolarRadiation_i.Resize(GridST2D);
        FileSolarRadiation_f.Resize(GridST2D);

        SolarRadiation_i.SetZero();
        SolarRadiation_f.SetZero();
        FileSolarRadiation_i.SetZero();
        FileSolarRadiation_f.SetZero();

	//----------------------------------

	/*** Background concentrations ***/
	Background_i.Resize(GridS2D, GridST2D);
	Background_f.Resize(GridS2D, GridST2D);
	FileBackground_i.Resize(GridS2D, GridST2D);
	FileBackground_f.Resize(GridS2D, GridST2D);
	
	Background_i.SetZero();
	Background_f.SetZero();
	FileBackground_i.SetZero();
	FileBackground_f.SetZero();

      }

    /*** Emissions rates ***/
    Emission_i.Resize(GridS2D, GridST2D);
    Emission_f.Resize(GridS2D, GridST2D);
    FileEmission_i.Resize(GridS2D, GridST2D);
    FileEmission_f.Resize(GridS2D, GridST2D);
	
    Emission_i.SetZero();
    Emission_f.SetZero();
    FileEmission_i.SetZero();
    FileEmission_f.SetZero();
    //-----------------------------------------

    /*** Street concentration ***/
    StreetConcentration.Resize(GridS2D, GridST2D);
    StreetConcentration.SetZero();

    /*** Deposition ***/

    StreetDryDepositionVelocity.Resize(GridS2D, GridST2D);
    StreetDryDepositionVelocity.SetZero(); 
    WallDryDepositionVelocity.Resize(GridS2D, GridST2D);
    WallDryDepositionVelocity.SetZero(); 
 
    
    // //StreetSurfaceDepositedMass.Resize(GridS2D, GridST2D);
	
    // StreetDryDepositionFlux.Resize(GridS2D, GridST2D);
    // WallDryDepositionFlux.Resize(GridS2D, GridST2D);
    // StreetDryDepositionRate.Resize(GridS2D, GridST2D); // ug/s
    // WallDryDepositionRate.Resize(GridS2D, GridST2D); // ug/s
	
    // //StreetSurfaceDepositedMass.SetZero();
    // StreetDryDepositionFlux.SetZero();
    // WallDryDepositionRate.SetZero();
    // WallDryDepositionFlux.SetZero();
    // StreetDryDepositionRate.SetZero();

    // /*** Scavenging ***/
    // ScavengingCoefficient.Resize(GridS2D, GridST2D);
    // StreetScavengingFlux.Resize(GridS2D, GridST2D);
    // StreetScavengingRate.Resize(GridS2D, GridST2D);

    // ScavengingCoefficient.SetZero();
    // StreetScavengingFlux.SetZero();
    // StreetScavengingRate.SetZero();

  }

  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template<class T>
  void StreetNetworkTransport<T>::Init()
  {
    this->SetCurrentDate(this->Date_min);
    InitStreet();
    InitIntersection();
    Allocate();
    ComputeStreetAngle();
    if (this->option_process["with_tree_aerodynamic"] or this->option_process["with_tree_deposition"])
      InitTree();
    if (this->option_process["with_tree_deposition"])
      ReadTreeParamDeposition();
  }

  //! Returns the name of a species with deposition velocities.
  /*!
    \param s species index in deposition velocities.
    \return The species name.
  */
  template<class T>
  string StreetNetworkTransport<T>::DepositionVelocityName(int s) const
  {
    return species_list_dep.at(s);
  }

  //! Returns the index of a species with deposition velocities.
  /*!
    \param s species index in deposition velocities.
    \return The index of a species in the species list.
  */
  template<class T>
  int StreetNetworkTransport<T>::DepositionVelocityGlobalIndex(int s) const
  {
   return this->GetSpeciesIndex(DepositionVelocityName(s));
  }
  
  //! Input data initialization.
  template<class T>
  void StreetNetworkTransport<T>::InitAllData()
  {

    string filename;
    /*** Meteo for the streets ***/
    if (this->option_process["with_local_data"])
      {

        if (this->option_process["with_transport"])
          {
            // *** Meteo for the streets ***/
            this->InitData("meteo",
                           "WindDirection",
                           FileWindDirection_i,
                           FileWindDirection_f,
                           this->current_date,
                           WindDirection_f);

            this->InitData("meteo",
                           "WindSpeed",
                           FileWindSpeed_i,
                           FileWindSpeed_f,
                           this->current_date,
                           WindSpeed_f);

            this->InitData("meteo",
                           "PBLH",
                           FilePBLH_i,
                           FilePBLH_f,
                           this->current_date,
                           PBLH_f);
        
            this->InitData("meteo",
                           "UST",
                           FileUST_i,
                           FileUST_f,
                           this->current_date,
                           UST_f);

            this->InitData("meteo",
                           "LMO",
                           FileLMO_i,
                           FileLMO_f,
                           this->current_date,
                           LMO_f);
          }

        
        this->InitData("meteo",
                       "SurfaceTemperature",
                       FileTemperature_i,
                       FileTemperature_f,
                       this->current_date,
                       Temperature_f);
        
        this->InitData("meteo",
                       "SurfacePressure",
                       FilePressure_i,
                       FilePressure_f,
                       this->current_date,
                       Pressure_f);

        this->InitData("meteo",
                       "Rain",
                       FileRain_i,
                       FileRain_f,
                       this->current_date,
                       Rain_f);

        this->InitData("meteo",
                       "SpecificHumidity",
                       FileSpecificHumidity_i,
                       FileSpecificHumidity_f,
                       this->current_date,
                       SpecificHumidity_f);

        // New meteo for gas deposition on tree leaves
        if (this->option_process["with_tree_deposition"])
          {
            this->InitData("meteo",
                           "SolarRadiation",
                           FileSolarRadiation_i,
                           FileSolarRadiation_f,
                           this->current_date,
                           SolarRadiation_f);
          }

        //---------------------------------------

        //  *** Meteo for the intersections ***/

        if (this->option_process["with_transport"])
          {

            this->InitData("meteo",
                           "WindDirectionInter",
                           FileWindDirectionInter_i,
                           FileWindDirectionInter_f,
                           this->current_date,
                           WindDirectionInter_f);
	
            this->InitData("meteo",
                           "WindSpeedInter",
                           FileWindSpeedInter_i,
                           FileWindSpeedInter_f,
                           this->current_date,
                           WindSpeedInter_f);
            
            this->InitData("meteo",
                           "PBLHInter",
                           FilePBLHInter_i,
                           FilePBLHInter_f,
                           this->current_date,
                           PBLHInter_f);
            
            this->InitData("meteo",
                           "LMOInter",
                           FileLMOInter_i,
                           FileLMOInter_f,
                           this->current_date,
                           LMOInter_f);
	
            this->InitData("meteo",
                           "USTInter",
                           FileUSTInter_i,
                           FileUSTInter_f,
                           this->current_date,
                           USTInter_f);
          }
      }
    
    //  *** Emissions ***/
    for (int i = 0; i < Ns_emis; i++)
      this->InitData("emission",
		     species_list_emis[i],
		     FileEmission_i,
		     FileEmission_f,
		     this->current_date,
		     i,
		     Emission_f);
    
    //! Background concentrations files.
    if (this->option_process["with_local_data"])
      {
	for (int i = 0; i < Ns_background; i++)
	  this->InitData("background_concentration",
			 species_list_background[i],
			 FileBackground_i,
			 FileBackground_f,
			 this->current_date,
			 i,
			 Background_f);
      }


    if (this->option_process["with_initial_condition"])
      {

        for (int i = 0; i < Ns_ic; i++)
          {
            string filename
              = this->input_files["initial_condition"](species_list_ic[i]);
            int index = this->GetSpeciesIndex(species_list_ic[i]);
            Data<T, 1>
              Concentration_tmp(&(StreetConcentration)(index, 0),
                                shape(total_nstreet));
            if (is_num(filename))
              Concentration_tmp.Fill(to_num<T>(filename));
            else
              FormatBinary<float>().Read(filename, Concentration_tmp);
          }
      }
  }


  //! Update the concentration arrays for each street.
  template<class T>
  void StreetNetworkTransport<T>::InitStreetConc()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        for (int s = 0; s < this->Ns; ++s)
          street->SetStreetConcentration(StreetConcentration(s, ist), s);
        
        ++ist;
      }
  }  


  
  //! Streets initialization.
  /*!
   */
  template<class T>
  void StreetNetworkTransport<T>::ReadStreetData()
  {
    /*** Street data ***/

    string line;
    vector<string> v;

    //! Get the input file name for the street data.
    this->config.SetSection("[street]");
    this->config.PeekValue("Street", file_street);

    total_nstreet = 0;
    ExtStream StreetStream(file_street);
    if (!StreetStream.is_open())
      throw string("File ") + file_street + " doesn't exist."; 

    //! Get the number of the streets.
    while (has_element(StreetStream))
      {
        StreetStream.GetLine(line);
        ++total_nstreet;
      }
    
    StreetStream.Rewind();
    //! Get the street data.
    id_street.resize(total_nstreet);
    begin_inter.resize(total_nstreet);
    end_inter.resize(total_nstreet);
    length.resize(total_nstreet);
    width.resize(total_nstreet);
    height.resize(total_nstreet);
    typo.resize(total_nstreet);
    typo = 0;
    Mean_length = 0.0;
    Mean_width = 0.0;
    Mean_height = 0.0;
    for (int i = 0; i < total_nstreet; ++i)
      {
        StreetStream.GetLine(line);
        v = split(line, ";");

        //! Check the number of values.
        if (v.size() != 7)
          throw string("Error: wrong number of values in street.dat ")
            + to_str(v.size() + " " + line); 
        
	id_street(i) = to_num<T>(v[0]);
        begin_inter(i) = to_num<T>(v[1]);
        end_inter(i) = to_num<T>(v[2]);
        length(i) = to_num<T>(v[3]);
        width(i) = to_num<T>(v[4]);
        height(i) = to_num<T>(v[5]);
	typo(i) = to_num<int>(v[6]);
	  	
        //! Check if a zero value exists.
        if (length(i) == 0.0)
          throw string("Street length is zero for the street ") + to_str(i + 1) + ".";
        if (width(i) == 0.0)
          throw string("Street width is zero for the street ") + to_str(i + 1) + ".";
        if (height(i) == 0.0)
          throw string("Builiding height is zero for the street ") + to_str(i + 1) + ".";

        // Calculation of the mean street length, width and bulding height for the whole street network
        Mean_length += length(i);
        Mean_width += width(i);
        Mean_height += height(i);
      }

    Mean_length /= total_nstreet; // Mean street length over the city
    Mean_width /= total_nstreet; // Mean street width over the city
    Mean_height /= total_nstreet; // Mean street height over the city
   
    // Calculation of the street length, width and building height to check street network homogeneity
    T variance_L = 0.0, variance_W = 0.0, variance_H = 0.0;
    T stdDeviation_L, stdDeviation_W, stdDeviation_H;
    for (int i = 0; i < total_nstreet; ++i)
      {
        variance_L += pow(length(i) - Mean_length,2);
        variance_W += pow(width(i) - Mean_width,2);
        variance_H += pow(height(i) - Mean_height,2);
      }
    variance_L = variance_L/total_nstreet;
    variance_W = variance_W/total_nstreet;
    variance_H = variance_H/total_nstreet;
    stdDeviation_L = sqrt(variance_L);
    stdDeviation_W = sqrt(variance_W);
    stdDeviation_H = sqrt(variance_H);
  }

  template<class T>
  void StreetNetworkTransport<T>::ReadTreeData()
  {
    /*** Tree data ***/

    string line;
    vector<string> v;

    //! Get the input file name for the tree data.
    this->config.SetSection("[street]");
    this->config.PeekValue("Tree", file_tree);

    total_ntree_street = 0;
    ExtStream TreeStream(file_tree);
    if (!TreeStream.is_open())
      throw string("File ") + file_tree + " doesn't exist.";

    //! Get the number of the streets in tree.dat.
    while (has_element(TreeStream))
      {
        TreeStream.GetLine(line);
        ++total_ntree_street;
      }
    if (total_nstreet != total_ntree_street)
      throw string("Missing streets in tree.dat : ") + to_str(total_nstreet) + " streets in street.dat vs " + to_str(total_ntree_street) + " in tree.dat.";
    TreeStream.Rewind();
    //! Get the tree data.
    id_street_tree.resize(total_nstreet);
    tree_height.resize(total_nstreet);
    trunk_height.resize(total_nstreet);
    tree_LAI.resize(total_nstreet);
    for (int i = 0; i < total_nstreet; ++i)
      {
        TreeStream.GetLine(line);
        v = split(line, ";");
        id_street_tree(i) = to_num<T>(v[0]);
        tree_height(i) = to_num<T>(v[1]);
        trunk_height(i) = to_num<T>(v[2]);
        tree_LAI(i) = to_num<T>(v[3]);
     }
  }

  //! Init tree data
  template<class T>
  void StreetNetworkTransport<T>::InitTree()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        street->SetTreeHeight(tree_height(ist));
        street->SetTrunkHeight(trunk_height(ist));
        street->SetTreeLAI(tree_LAI(ist));
        if (tree_height(ist) > street->GetHeight())
          cout << "Tree crown height exceeds the building height in street number " << id_street_tree(ist) << "." << endl;
        ++ist;
      }
  }

  // Read the parameter file to get model parameters for dry deposition on tree leaves
  template<class T>
  void StreetNetworkTransport<T>::ReadTreeParamDeposition()
  {
    string line;
    vector<string> v;

    ExtStream NlucTreeDep(file_tree_deposition);
    NlucTreeDep.GetLine(line);
    v = split(line, "\t ");
    Nluc = v.size();

    Data<T, 1> LUC(Nluc);
    Data<T, 1> RcutSO2_d(Nluc);
    Data<T, 1> RcutO3_d(Nluc);
    Data<T, 1> Rsmin(Nluc);
    Data<T, 1> Radius(Nluc);
    Data<T, 1> Alpha(Nluc);
    Data<T, 1> Gamma(Nluc);

    ifstream ParamTreeDep;
    ParamTreeDep.open(file_tree_deposition.c_str());
    FormatText InputParam;

    InputParam.Read(ParamTreeDep, LUC);
    InputParam.Read(ParamTreeDep, RcutSO2_d);
    InputParam.Read(ParamTreeDep, RcutO3_d);
    InputParam.Read(ParamTreeDep, Rsmin);
    InputParam.Read(ParamTreeDep, Radius);
    InputParam.Read(ParamTreeDep, Alpha);
    InputParam.Read(ParamTreeDep, Gamma);

    ParamTreeDep.close();

    for (int i = 0; i < Nluc; i++)
      {
        if (tree_LUC == int(LUC(i)))
          {
            // Get the model parameters for gas deposition on tree leaves
            RcutSO2_d0 = RcutSO2_d(i);
            RcutO3_d0 = RcutO3_d(i);
            Rsmin_tree = Rsmin(i);
            // Get the Zhang model parameters Radius, Alpha and Gamma for aerosol deposition on tree leaves
            Radius_tree = Radius(i);
            alpha_tree = Alpha(i);
            gamma_tree = Gamma(i);
          }
      }
 }

  //! Streets initialization.
  /*!
   */
  template<class T>
  void StreetNetworkTransport<T>::InitStreet()
  {

    Array<T, 1> init_conc(this->Ns);
    init_conc = 0.0;
    for (int i = 0; i < total_nstreet; ++i)
      {
        Street<T>* street = 
          new Street<T>(id_street(i),
			begin_inter(i),
			end_inter(i),
                        length(i),
			width(i),
			height(i),
			typo(i),
                        this->Ns);
        for (int s = 0; s < this->Ns; s++)
            street->SetStreetConcentration(init_conc(s), s); 

        StreetVector.push_back(street);
      }
    current_street = StreetVector.begin();
  }
  
  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T>
  void StreetNetworkTransport<T>::InitStep()
  {
    BaseModel<T>::InitStep();
    
    //  *** Emissions ***/
    for (int i = 0; i < Ns_emis; i++)
      this->UpdateData("emission",
		       species_list_emis[i],
		       FileEmission_i,
		       FileEmission_f,
		       i,
		       Emission_f);

    
    if (this->option_process["with_local_data"])
      {
	//! Background concentrations files.
	for (int i = 0; i < Ns_background; i++)
	  this->UpdateData("background_concentration",
			   species_list_background[i],
			   FileBackground_i,
			   FileBackground_f,
			   i,
			   Background_f);


        //! Meteo for streets.
        if (this->option_process["with_transport"])
          {
            this->UpdateData("meteo",
                             "WindDirection",
                             FileWindDirection_i,
                             FileWindDirection_f,
                             WindDirection_f);

            this->UpdateData("meteo",
                             "WindSpeed",
                             FileWindSpeed_i,
                             FileWindSpeed_f,
                             WindSpeed_f);
	
            this->UpdateData("meteo",
                             "PBLH",
                             FilePBLH_i,
                             FilePBLH_f,
                             PBLH_f);
	
            this->UpdateData("meteo",
                             "UST",
                             FileUST_i,
                             FileUST_f,
                             UST_f);
	
            this->UpdateData("meteo",
                             "LMO",
                             FileLMO_i,
                             FileLMO_f,
                             LMO_f);
          }

        
	this->UpdateData("meteo",
			 "SurfaceTemperature",
			 FileTemperature_i,
			 FileTemperature_f,
			 Temperature_f);
	
	this->UpdateData("meteo",
			 "SurfacePressure",
			 FilePressure_i,
			 FilePressure_f,
			 Pressure_f);
	
	this->UpdateData("meteo",
			 "Rain",
			 FileRain_i,
			 FileRain_f,
			 Rain_f);

    
        this->UpdateData("meteo",
                         "SpecificHumidity",
                         FileSpecificHumidity_i,
                         FileSpecificHumidity_f,
                         SpecificHumidity_f);

        if (this->option_process["with_tree_deposition"])
          {
            this->UpdateData("meteo",
                             "SolarRadiation",
                             FileSolarRadiation_i,
                             FileSolarRadiation_f,
                             SolarRadiation_f);
          }
	
     	//--------------------------------------------

        //!  Meteo for intersections.
        if (this->option_process["with_transport"])
          {

            this->UpdateData("meteo",
                             "WindDirectionInter",
                             FileWindDirectionInter_i,
                             FileWindDirectionInter_f,
                             WindDirectionInter_f);
	
            this->UpdateData("meteo",
                             "WindSpeedInter",
                             FileWindSpeedInter_i,
                             FileWindSpeedInter_f,
                             WindSpeedInter_f);
	
            this->UpdateData("meteo",
                             "PBLHInter",
                             FilePBLHInter_i,
                             FilePBLHInter_f,
                             PBLHInter_f);
	
            this->UpdateData("meteo",
                             "USTInter",
                             FileUSTInter_i,
                             FileUSTInter_f,
                             USTInter_f);
	
            this->UpdateData("meteo",
                             "LMOInter",
                             FileLMOInter_i,
                             FileLMOInter_f,
                             LMOInter_f);
          }
      }

    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	Array<T,1> emission_rate(this->Ns);
        //! Set the emission data
        emission_rate = 0.0; 
        for (int i = 0; i < this->Ns; ++i)
          for (int j = 0; j < Ns_emis; ++j)
            if (this->species_list[i] == species_list_emis[j])
	      {
		if(Emission_f(j,st)*1 != Emission_f(j,st))
		  throw string("NaN in emission data");
		emission_rate(i) = Emission_f(j, st);
	      }
        street->SetEmission(emission_rate);

        if (this->option_process["with_local_data"])
          {
            //! Set the meteo data
	    street->SetMeteoTransport(WindDirection_f(st),
	    			      WindSpeed_f(st),
	    			      PBLH_f(st),
	    			      UST_f(st),
	    			      LMO_f(st),
	    			      Pressure_f(st),
	    			      Temperature_f(st),
	    			      Rain_f(st),
				      SpecificHumidity_f(st),
                                      // New meteo for gas deposition on tree leaves
                                      SolarRadiation_f(st));

	    //relative humidity
	    T relative_humidity_ = 0.0;
	    _compute_relative_humidity(&SpecificHumidity_f(st),
				       &Temperature_f(st),
				       &Pressure_f(st),
				       &relative_humidity_);
	    if(relative_humidity_ *1. != relative_humidity_)
	      throw string("NaN in relative_humidity in street " +to_str(st) + " : " + to_str(relative_humidity_) );
	    relative_humidity_ =
	      min(max(relative_humidity_, 0.0), 0.95);
	    street->SetRelativeHumidity(relative_humidity_);

            //! Set the background concentration data
            for (int i = 0; i < this->Ns; ++i)
              for (int j = 0; j < Ns_background; ++j)
                if (this->species_list[i] == species_list_background[j])
		  {
		    if(Background_f(j, st) * 1 != Background_f(j, st))
		      throw string("NaN in background concentration");
		    street->SetBackgroundConcentration(Background_f(j,st), i);
            if (!this->option_process["with_initial_condition"])
              {
                //LL: Concentration at streets equal to Cbg in first time step
                if(this->current_date == this->Date_min)
                  street->SetStreetConcentration(Background_f(j,st), i);
              }

		  }
          }
        ++st; 
      }

    //! Initilialize the inflow rate.
    InitInflowRate();

    //! Initialize the flux exchange matrix at the intersections.
    int inter = 0;
    for (typename vector<Intersection<T>* >::iterator iter = IntersectionVector.begin(); iter != IntersectionVector.end(); iter++)
      {
        Intersection<T>* intersection = *iter;
        int nstreet_inter = intersection->GetNStreet();
        Array<T, 2> zeromatrix(nstreet_inter + 1, nstreet_inter + 1);
        zeromatrix = 0.0;
        intersection->SetFluxMatrix(zeromatrix);
        intersection->SetGaussianMatrix(zeromatrix);

        //! Set the meteo data
        if (this->option_process["with_local_data"])
          {
            intersection->SetMeteo(WindDirectionInter_f(inter),
                                   WindSpeedInter_f(inter),
                                   PBLHInter_f(inter),
                                   USTInter_f(inter),
                                   LMOInter_f(inter));
          }
        ++inter;
      }
  }

  //! Sets the current street index and iterator to a given street.
  /*! It sets current street index to the index given as parameter, and
    current street points to the corresponding street in StreetVector. It returns
    the current street.
    \param index street index.
  */
  template<class T>
  void StreetNetworkTransport<T>::SetCurrentStreet(int index)
  {
    int street_index = 0;
    current_street = StreetVector.begin();
    if (street_index < index)
      while (street_index != index)
	{
	  if (current_street != StreetVector.end())
	    {
	      street_index ++;
	      current_street ++;
	    }
	  else
	    throw string("Street index: ") + to_str<int>(index)
	      + " out of range.";
	}
    else if (street_index > index)
      {
	street_index = 0;
	current_street = StreetVector.begin();
	while (street_index != index)
	  {
	    if (current_street != StreetVector.end())
	      {
		street_index ++;
		current_street ++;
	      }
	    else
	      throw string("Street index: ") + to_str<int>(index)
		+ " out of range.";
	  }
      }
  }

  //! Sets the current intersection index and iterator to a given intersection.
  /*! It sets current intersection index to the index given as parameter, and
    current intersection points to the corresponding intersection in IntersectionVector. It returns
    the current intersection.
    \param index intersection index.
  */
  template<class T>
  void StreetNetworkTransport<T>::SetCurrentIntersection(int index)
  {
    int intersection_index = 0;
    current_intersection = IntersectionVector.begin();
    while (intersection_index != index)
      {
        if (current_intersection != IntersectionVector.end())
          {
            intersection_index ++;
            current_intersection ++;
          }
        else
          throw string("Intersection index: ") + to_str<int>(index)
            + " out of range.";
      }
  }

  //! Returns the quantity of species s in a given puff.
  /*!
    \return The quantity of species s in a given puff.
  */
  template<class T>
  T StreetNetworkTransport<T>::GetStreetQuantity(int index, int s)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;
    return (*current_street)->GetStreetQuantity(s);
  }

  //! Returns the emission rate of species s in a given street index.
  /*!
    \return The emission rate of species s in a given street index.
  */
  template<class T>
  T StreetNetworkTransport<T>::GetStreetEmissionRate(int index, int s)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;
    return street->GetEmission(s);
  }

  //! Returns the street height.
  /*!
    \return The street height.
  */
  template<class T>
  T StreetNetworkTransport<T>::GetStreetHeight(int index)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;
    return street->GetHeight();
  }

  //! Returns the street length.
  /*!
    \return The street length.
  */
  template<class T>
  T StreetNetworkTransport<T>::GetStreetLength(int index)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;
    return street->GetLength();
  }

  //! Returns the street ID.
  /*!
    \return The street ID.
  */
  template<class T>
  int StreetNetworkTransport<T>::GetStreetID(int index)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;
    return street->GetStreetID();
  }


  //! Returns the street canopy volume.
  /*!
    \return The street canopy volume.
  */
  template<class T>
  T StreetNetworkTransport<T>::GetStreetVolume(int index)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;
    return (*current_street)->GetStreetVolume();
  }


  //! Returns the mass transfer rate to the background.
  /*!
    \return The mass transfer rate to the background.
  */
  template<class T>
  T StreetNetworkTransport<T>::GetMassTransferBackground(int index, int s)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;

    //! positive value: mass transfer from street to background
    //! negative value: mass transfer from background to street
    T quantity_delta = street->GetStreetQuantityDelta(s);

    T mass_transfer_to_background;
    
    T roof_to_background, intersection_to_background, intersection_from_background;
    roof_to_background = street->GetMassfluxRoofToBackground(s);
    
    intersection_to_background = street->GetMassfluxToBackground(s);
    intersection_from_background = street->GetMassfluxFromBackground(s);

    mass_transfer_to_background = roof_to_background + intersection_to_background
      - intersection_from_background; // ug/s

    return mass_transfer_to_background;
  }


  //! Returns the mass flux balance.
  /*!
    \return The mass flux balance.
  */
  template<class T>
  T StreetNetworkTransport<T>::GetMassFluxExchange(int index, int s)
  {
    SetCurrentStreet(index);
    Street<T>* street = *current_street;
    
    T temp1, temp3, temp4;
    temp1 = street->GetMassfluxRoofToBackground(s);
    temp3 = street->GetMassfluxToBackground(s);
    temp4 = street->GetMassfluxFromBackground(s);

    return temp1 + temp3 - temp4;
  }


  //! Erases a street from the street list.
  template<class T>
  void StreetNetworkTransport<T>::EraseStreet()
  {
    delete *current_street;
    StreetVector.erase(current_street);
  }

  //! Clears the street list.
  template<class T>
  void StreetNetworkTransport<T>::ClearStreetVector()
  {
    for (typename vector<Street<T>* >:: iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      delete (*iter);
    StreetVector.clear();
  }

  //! Intersections initialization.
  template<class T>
  void StreetNetworkTransport<T>::InitIntersection()
  {
    /*** Intersection data ***/

    //! Get the input file for the intersection data
    this->config.SetSection("[street]");
    this->config.PeekValue("Intersection", file_intersection);

    string line;
    vector<string> v;
    nintersection = 0;

    ExtStream IntersectionStream(file_intersection);
    if (!IntersectionStream.is_open())
      throw string("File ") + file_intersection + " doesn't exist."; 

    //! Get the number of intersections.
    while (has_element(IntersectionStream))
      {
        IntersectionStream.GetLine(line);
        ++nintersection;
      }
    IntersectionStream.Rewind();

    //! Read the intersection data.
    id_inter.resize(nintersection);
    x_inter.resize(nintersection);
    y_inter.resize(nintersection);
    nstreet.resize(nintersection);
    is_virtual.resize(nintersection);
    maxnstreet = 10;
    street_list.resize(nintersection, maxnstreet);
    street_list = 0;
    for (int i = 0; i < nintersection; ++i)
      {
        IntersectionStream.GetLine(line);
        v = split(line, ";");
        id_inter(i) = to_num<T>(v[0]);
        x_inter(i) = to_num<T>(v[1]);
        y_inter(i) = to_num<T>(v[2]);
        nstreet(i) = to_num<T>(v[3]);
        int ndata = v.size();
        if (ndata != nstreet(i) + 4)
          throw string("Number of data is invalid in ") + file_intersection + 
            + ". " + to_str(nstreet(i) + 4) + " are needed but " 
            + to_str(ndata) + " are given. Intersection " + to_str(id_inter(i)) ;

        for (int j = 0; j < nstreet(i); ++j)
          street_list(i, j) = to_num<T>(v[4 + j]);
        if (nstreet(i) == 1)
          is_virtual(i) = true;
        else
          is_virtual(i) = false;
      }


    for (int i = 0; i < nintersection; ++i)
      {
        Array<int, 1> street_list_inter(nstreet(i));
        for (int j = 0; j < nstreet(i); j++)
          street_list_inter(j) = street_list(i, j);
        Array<T, 2> flux_matrix(nstreet(i) + 1, nstreet(i) + 1);
        flux_matrix = 0.0;
        Array<T, 2> gaussian_matrix(nstreet(i) + 1, nstreet(i) + 1);
        gaussian_matrix = 0.0;
        Intersection<T>* intersection = 
          new Intersection<T>(id_inter(i), x_inter(i), y_inter(i),
                              nstreet(i), street_list_inter, is_virtual(i),
                              flux_matrix, gaussian_matrix);
        IntersectionVector.push_back(intersection);
      }
    current_intersection = IntersectionVector.begin();
  }

  //! Erases a intersection from the intersection list.
  template<class T>
  void StreetNetworkTransport<T>::EraseIntersection()
  {
    delete *current_intersection;
    IntersectionVector.erase(current_intersection);
  }

  //! Clears the intersection list.
  template<class T>
  void StreetNetworkTransport<T>::ClearIntersectionVector()
  {
    for (typename vector<Intersection<T>* >:: iterator iter = IntersectionVector.begin();
         iter != IntersectionVector.end(); iter++)
      delete (*iter);
    IntersectionVector.clear();
  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T>
  void StreetNetworkTransport<T>::CheckConfiguration()
  {
    if(sub_delta_t_min > this->Delta_t)
      throw string("Error! sub_delta_t_min can not be higher than simulation time step.");
    
    /*** Meteo ***/
    if (this->option_process["with_local_data"])
      {

      	if(this->option_process["with_transport"])
          {
            //! Check the meteo data for the streets.
            if (this->input_files["meteo"]("WindDirection").empty())
              throw "WindDirection is needed but no input data file was provided.";
            if (this->input_files["meteo"]("WindSpeed").empty())
              throw "WindSpeed is needed but no input data file was provided.";
            if (this->input_files["meteo"]("PBLH").empty())
              throw "PBLH is needed but no input data file was provided.";
            if (this->input_files["meteo"]("UST").empty())
              throw "UST is needed but no input data file was provided.";
            if (this->input_files["meteo"]("LMO").empty())
              throw "LMO is needed but no input data file was provided.";

            //! Check the meteo data for the intersections.
            if (this->input_files["meteo"]("WindDirectionInter").empty())
              throw "WindDirectionInter is needed but no input data file was provided.";
            if (this->input_files["meteo"]("WindSpeedInter").empty())
              throw "WindSpeedInter is needed but no input data file was provided.";
            if (this->input_files["meteo"]("PBLHInter").empty())
              throw "PBLHInter is needed but no input data file was provided.";
            if (this->input_files["meteo"]("USTInter").empty())
              throw "USTInter is needed but no input data file was provided.";
            if (this->input_files["meteo"]("LMOInter").empty())
              throw "LMOInter is needed but no input data file was provided.";
          }

        //! Data for deposition 
        if(this->option_process["with_deposition"])
          {
            if (this->input_files["meteo"]("SurfaceTemperature").empty())
              throw "Temperature is needed but no input data file was provided.";
            if (this->input_files["meteo"]("SurfacePressure").empty())
              throw "Pressure is needed but no input data file was provided.";
          }
    
        if(this->option_process["with_scavenging"] and
           this->scavenging_model != "none")
          {
            if (this->input_files["meteo"]("SurfaceTemperature").empty())
              throw "Temperature is needed but no input data file was provided.";
            if (this->input_files["meteo"]("SurfacePressure").empty())
              throw "Pressure is needed but no input data file was provided.";
            if (this->input_files["meteo"]("Rain").empty())
              throw "Rain is needed but no input data file was provided.";
          }
      }
    
    if (this->option_process["with_scavenging"] and
        this->scavenging_model == "none")
      throw "Warning!! Option With_scavenging = yes and scavenging_model = none";   
    if (!this->option_process["with_deposition"]
        && this->option_process["collect_dry_flux"])
      throw string("Dry deposition fluxes cannot be collected") +
        " without deposition.";

    if (this->option_process["collect_wet_flux"]
        && (!this->option_process["with_scavenging"]))
      throw string("Wet deposition fluxes cannot be collected") +
        " without scavenging.";
    
  }
  

  //! Performs one step forward.
  template<class T>
  void StreetNetworkTransport<T>::Forward()
  {

    InitMeteo();

    ComputeMassBalance();

    SetStreetConcentration();
    SetStreetDryDepositionVelocity();
    SetWallDryDepositionVelocity();
    this->AddTime(this->Delta_t);
    this->step++;
  }

  //! Calls the functions for the transport.
  template<class T>
  void StreetNetworkTransport<T>::InitMeteo()    
  {
    is_stationary = false;

    //! Compute the wind speed in the street-canyon.    
    ComputeUstreet();
    ComputeSigmaW();
    ComputeTransferVelocity();
    ComputeWindDirectionFluctuation();

    //! Compute deposition velocity
    if(this->option_process["with_deposition"])
      ComputeDryDepositionVelocities();

    //! Compute scavenging coefficient
    if(this->option_process["with_scavenging"])
      ComputeScavengingCoefficient();
    
    SetInitialStreetConcentration();

  }


  //! Compute mass balance.
  template<class T>
  void StreetNetworkTransport<T>::ComputeMassBalance()    
  {
    
    if (this->option_process["with_stationary_hypothesis"])
      {
	int niter = 0;
	const int niter_max = nintersection;
	while (niter < niter_max and (!is_stationary))
	  {
	    InitInflowRate();
	    // cout << " ----> MUNICH Iteration No " << niter << endl;
	    ComputeInflowRateExtended();
	    //! Compute the concentrations in the street-canyon.
	    ComputeStreetConcentration();
	    IsStationary(is_stationary);
	    ++niter;
	  }
	if (!is_stationary)
	  throw("Error: stationarity is not achieved. Please increase the number of iterations.");
      }
    else //non stationary calcul
      {
	InitInflowRate();
	ComputeInflowRateExtended();
	ComputeStreetConcentrationNoStationary();
      }
  }
  
  
  //Set new street concentration and deposition
  template<class T>
  void StreetNetworkTransport<T>::SetStreetOutput()
  {
    SetStreetConcentration();
    if(this->option_process["with_deposition"])
      {
        SetStreetDryDepositionVelocity();
        SetWallDryDepositionVelocity();
        //SetStreetDryDepositionFlux();
        //SetWallDryDepositionFlux();
        //SetStreetDryDepositionRate(); // ug/s
        //SetWallDryDepositionRate(); // ug/s	
      }
    if(this->option_process["with_scavenging"])
      {
        SetStreetScavengingFlux();
        SetStreetScavengingRate();
      }
  }

  //! Checks if the stationarity is acheved for every streets.
  template<class T>
  void StreetNetworkTransport<T>::IsStationary(bool& is_stationary)
  {
    is_stationary = true;
    //! Get the mass flux balance for the atmosphere.
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        if (!street->GetStationary())
          is_stationary = false;
      }
  }


  //! Compute fluxes at the intersection.
  /*!
    \param wind_dir_inter wind direction at the intersection
    \param intersection intersection object
  */
  template<class T>
  void StreetNetworkTransport<T>::ComputeIntersection(T wind_dir_inter, 
                                                      Intersection<T>* intersection)
  {
    //! Streets which are connected to the intersection.
    int nstreet_inter = intersection->GetNStreet(); 
    Array<T, 2> extended_matrix(nstreet_inter + 1, nstreet_inter + 1);
    extended_matrix = 0.0;
    Array<int, 1> street_list_inter(nstreet_inter);
    street_list_inter = intersection->GetStreetList();

    if (!intersection->IsVirtual())
      {
        for (int j = 0; j < nstreet_inter; ++j)
          {
            int end_inter_id = 0;
            int begin_inter_id = 0;
            for (typename vector<Street<T>* >::iterator iter2 = StreetVector.begin(); 
                 iter2 != StreetVector.end(); iter2++)
              {
                Street<T>* street = *iter2;
                if (street_list_inter(j) == street->GetStreetID())
                  {
                    StreetVectorInter.push_back(street);
                    break;
                  }
              }
            
            /* When the street angle was calculated, the begining intersection of 
               the street was considered as the origin point.
               Therefore if the intersection is the ending intersection of 
               the street the street angle should be corrected by 
               addtion or subtraction by pi.
            */
            Street<T>* street = StreetVectorInter[j];
            if (intersection->GetID() == street->GetEndIntersectionID())
              {
                if (street->GetStreetAngle() >= pi)
                  street->SetStreetAngleIntersection(street->GetStreetAngle() - pi);
                else
                  street->SetStreetAngleIntersection(street->GetStreetAngle() + pi);
              }
            else if (intersection->GetID() == 
                     street->GetBeginIntersectionID())
              street->SetStreetAngleIntersection(street->GetStreetAngle());
            else
              throw string("Error: unknown ID") + to_str(intersection->GetID());
          } // nstreet_inter

        ComputeIntersectionFlux(extended_matrix, wind_dir_inter);

        intersection->SetFluxMatrix(extended_matrix);

        //! Clear the street vector for the intersection. 
        StreetVectorInter.clear();
      } // IsVirtual: false
    else // End-node case
      {
        // Input variables
        if (nstreet_inter != 1)
          throw string("Error: the virtual intersection corresponds to a single street: ") + to_str(nstreet_inter) + " given.";
        Array<int, 1> street_list_inter(nstreet_inter);
        street_list_inter = intersection->GetStreetList();
       
        // Virtual intersections
        for (typename vector<Street<T>* >::iterator iter2 = StreetVector.begin(); iter2 != StreetVector.end(); iter2++)
          {
            Street<T>* street = *iter2;
            if (street_list_inter(0) == street->GetStreetID())
              {
                if (intersection->GetID() == street->GetEndIntersectionID())
                  {
                    if (street->GetStreetAngle() >= pi)
                      street->SetStreetAngleIntersection(street->GetStreetAngle() - pi);
                    else
                      street->SetStreetAngleIntersection(street->GetStreetAngle() + pi);
                  }
                else if (intersection->GetID() == 
                         street->GetBeginIntersectionID())
                  street->SetStreetAngleIntersection(street->GetStreetAngle());
                else
                  throw string("Error: unknown ID") + to_str(intersection->GetID());

                T dangle = abs(street->GetStreetAngleIntersection() - wind_dir_inter);
                T flux = street->GetHeight() * 
                  street->GetWidth() * street->GetStreetWindSpeed();

                if (dangle > (pi / 2.0) and dangle < (pi * 3.0 / 2.0))
                  {
                    //! Incoming flow to the atmosphere, 
                    //! i.e., outgoing flow from the street. 
                    extended_matrix(1, 0) = flux;
                  }
                else
                  {
                    //! Outgoing flow from the atmosphere,
                    //! i.e., incoming flow to the street.
                    extended_matrix(0, 1) = flux;
                  }
                intersection->SetFluxMatrix(extended_matrix);
                break;
              }
          }
      } // End-node case

  }

  //! Checks whether a species is scavenged.
  /*!
    \param s species global index.
    \return True if the species is scavenged, false otherwise.
  */
  template<class T>
  bool StreetNetworkTransport<T>
  ::HasScavenging(int s) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(),
		this->GetSpeciesName(s)) != species_list_scav.end();
  }

  //! Checks whether a species is scavenged.
  /*!
    \param name species name.
    \return True if the species is scavenged, false otherwise.
  */
  template<class T>
  bool StreetNetworkTransport<T>
  ::HasScavenging(string name) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(), name)
      != species_list_scav.end();
  }
  
  //! (Re)Initialize the inflow rate.
  template<class T>
  void StreetNetworkTransport<T>::InitInflowRate()
  {
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); 
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
          {
            street->SetInflowRate(0.0, s);
            street->SetOutgoingFlux(0.0);
            street->SetMassfluxFromBackground(0.0, s);
            street->SetMassfluxToBackground(0.0, s);
	    street->SetIncomingFlux(0.0);
          }
      }
  }


//! Returns the index in scavenged species of a given species.
  /*!
    \param name species name.
    \return The species index in scavenged species.
  */
  template<class T>
  bool StreetNetworkTransport<T>
  ::ScavengingIndex(string name) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(), name)
      - species_list_scav.begin();
  }

  //! Returns the index in scavenged species of a given species.
  /*!
    \param s species global index.
    \return The species index in scavenged species.
  */
  template<class T>
  bool StreetNetworkTransport<T>
  ::ScavengingIndex(int s) const
  {
    return find(species_list_scav.begin(), species_list_scav.end(),
		this->GetSpeciesName(s)) - species_list_scav.begin();
  }

  // Sets infinite resistance values to 1e+30.
  template<class T>
  void StreetNetworkTransport<T>::Infinity(Data<T, 1>& Data_info)
  {
    int Nc = Data_info.GetLength(0);
    for(int c = 0; c < Nc; c++)
      if (Data_info(c) == 9999)
	Data_info(c) = 1e+30;
  }

  // Sets infinite resistance values to 1e+30.
  template<class T>
  void StreetNetworkTransport<T>::Cut(Data<T, 1>& Data_info)
  {
    int Nc = Data_info.GetLength(0);
    for(int c = 0; c < Nc; c++)
      {
	Data_info(c) = max(Data_info(c), 0.0);
	Data_info(c) = min(Data_info(c), 1.0e+30);
      }
  }

  // Sets infinite resistance values to 1e+30.
  template<class T>
  void StreetNetworkTransport<T>::Cut(T& Data_info)
  {
    Data_info = max(Data_info, 0.0);
    Data_info = min(Data_info, 1.0e+30);
  }

  
  //! Compute surface friction velocity for deposisiton velocity according to the wind speed options
  template<class T>
  T StreetNetworkTransport<T>::ComputeUstarAboveSurface(T z1, T H, T W, T z0s, T ustar_city)
  {
    T ustar_z1;
    T Uh = 0.0;
    T Um = 0.0;
    if (option_uH == "Sirane")
      {
        T delta = min(H, W / 2);
        T C = ComputeSiraneC(z0s, delta);
        Um = ComputeSiraneUm(C, ustar_city, z0s);
        if (option_ustreet == "Sirane")
          ustar_z1 = ComputeSiraneUstarProfile(z1, H, C, delta, Um);
        else if (option_ustreet == "Wang")
          {
            T sH = ComputeWangsH(H, W, 0.0, 0.0, 0.0);
            ustar_z1 = ComputeWangUstarProfile(z1, H, W, z0s, sH, ustar_city, 0.0, 0.0);
          }
        else if (option_ustreet == "Exponential")
          {
            ustar_z1 = ComputeExpUstarProfile(z1, H, W, ustar_city);
          }
        else
          throw("Wrong option given. Choose Sirane, Wang or Exponential.");
      }
    else if (option_uH == "Macdonald")
      {
        ComputeMacdonaldProfile();
        if (option_ustreet == "Sirane")
          {
            if  (H < (d_city + z0_city))
              Uh = 0.0;
            else
              Uh = ComputeMacdonaldUH(H, ustar_city, d_city, z0_city);
            T delta = min(H, W / 2);
            T C = ComputeSiraneC(z0s, delta);
            T Uh_Um_factor = ComputeUmtoUHFactor(C, z0s);
            Um = Uh / Uh_Um_factor;
            ustar_z1 = ComputeSiraneUstarProfile(z1, H, C, delta, Um);
          }
        else if (option_ustreet == "Wang")
          {
            T sH = ComputeWangsH(H, W, 0.0, 0.0, 0.0);
            ustar_z1 = ComputeWangUstarProfile(z1, H, W, z0s, sH, ustar_city, 0.0, 0.0);
          }
        else if (option_ustreet == "Exponential")
          ustar_z1 = ComputeExpUstarProfile(z1, H, W, ustar_city);
        else
          throw("Wrong option given. Choose Sirane, Wang or Exponential.");
      }
    else
      throw("Wrong option given. Choose Sirane or Macdonald.");
  return ustar_z1;
  }

  //! Compute dry deposition of gas species for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::ComputeDryDepositionVelocities()
  {
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        T H = street->GetHeight();
        T W = street->GetWidth();
        T temperature_ = street->GetTemperature();
        T temperature_celsius = temperature_ - 273.15;

        // Compute wall and ground surface friction velocity
        T z1 = z0s; // altitude above the surface to compute the friction velocty (m)
        T ustar_city = street->GetStreetUstar();
        T ustar_surface = ComputeUstarAboveSurface(z1, H, W, z0s, ustar_city);
 
        /*** Input data for Rb (Quasi-laminar sublayer resistance or diffusion resistance) ***/ 
        Array<T, 1> gas_phase_diffusivity_(Ns_dep), Sc(Ns_dep);
        Array<T, 1> alpha_(Ns_dep), beta_(Ns_dep);
        Array<T, 1> Rb(Ns_dep);

        for (int s = 0; s < Ns_dep; s++)
          {
            /*** Rb ***/ 
            // Quasi-laminar sublayer resistance (diffusion resistance)
            gas_phase_diffusivity_(s) = gas_phase_diffusivity[DepositionVelocityName(s)];
            if (gas_phase_diffusivity_(s) == 0.0)
                throw string("Error: diffusivity is required for ") +
                  DepositionVelocityName(s) + ".";            
            Sc(s) = ComputeSchmidtNumber(nu, gas_phase_diffusivity_(s));
            Rb(s) = ComputeBoundaryLayerResistance(ustar_surface, Sc(s));
 
            /*** Rg ***/
            // Surface resistance (ground resistance)
            alpha_(s) = alpha[DepositionVelocityName(s)];
            beta_(s) = beta[DepositionVelocityName(s)];
            T Rg_O3 = 500.; // values for dry surfaces
            T Rg_SO2 = 300.;
            T Rg;
            if (DepositionVelocityName(s) == "O3")
              Rg = Rg_O3;
            else if (DepositionVelocityName(s) == "SO2")
              Rg = Rg_SO2;
            else
              Rg = ComputeGroundResistance(Rg_SO2, Rg_O3, alpha_(s), beta_(s), temperature_celsius);

	    Cut(Rg);
            T Req = Rb(s) + Rg;
	    
            if (Req <= 0.0)
              throw string("Error in deposition, R_street: ")
                + to_str(Req);
            // Walls and street ground have the same roughness length and surface friction velocity
            // so the dry deposition velocity is equal on both surfaces.
            T street_dry_deposition_velocity = 1. / Req;
            street->SetStreetDryDepositionVelocity(street_dry_deposition_velocity, 
                                                   DepositionVelocityGlobalIndex(s));
            street->SetWallDryDepositionVelocity(street_dry_deposition_velocity, 
                                                 DepositionVelocityGlobalIndex(s));
          }

          // Compute gas dry deposition on tree leaves
        T tree_dry_deposition_velocity = 0.0;
        T LAI, hmax, htrunk = 0.0;
        if (this->option_process["with_tree_deposition"])
          {
            hmax = street->GetTreeHeight();
            htrunk = street->GetTrunkHeight();
            LAI = street->GetTreeLAI();
            if (LAI != 0.0 and hmax != 0.0 and htrunk != 0.0)
              {
                // Compute friction velocity in the tree crown
                T sH = ComputeWangsH(H, W, hmax, LAI, Cdt);
                T hm = htrunk + (hmax - htrunk) / 2; // height of the middle of the tree crown
                T ustar_htree = ComputeWangUstarProfile(hm, H, W, z0s, sH, ustar_city, LAI, Cdt);
                Array<T, 1> henry_constant_(Ns_dep), reactivity_(Ns_dep);

                for (int s = 0; s < Ns_dep; s++)
                  {
                    /*** Rb_tree ***/
                    // Quasi-laminar sublayer resistance (diffusion resistance)
                    T Rb_tree = ComputeBoundaryLayerResistance(ustar_htree, Sc(s));

                    /*** Rmes ***/
                    // Compute mesophyll resistance
                    henry_constant_(s) = henry_constant[DepositionVelocityName(s)];
                    reactivity_(s) = reactivity[DepositionVelocityName(s)];
                    T Rmes = ComputeMesophyllResistance(henry_constant_(s), reactivity_(s));
                    Cut(Rmes);

                    /*** Rcut ***/
                    // Compute cuticular resistance
                    T Rcut;
                    if (option_Rcut == "zhang")
                      {
                        T RH = street->GetRelativeHumidity();
                        Rcut = ComputeZhangCuticularResistance(RcutSO2_d0, RcutO3_d0, alpha_(s), beta_(s), temperature_celsius, RH, ustar_htree, LAI);
                      }
                    else if (option_Rcut == "wesely")
                      Rcut = ComputeWeselyCuticularResistance(henry_constant_(s), reactivity_(s), temperature_celsius, LAI);
                    else
                      throw("Wrong cuticular resistance option given. Choose 'zhang' or 'wesely'.");
                    Cut(Rcut);

                    /*** Rsto ***/
                    // Compute stomatal resistance
                    T DmH2O = 0.260; // Diffusion coefficient of water vapor in air (cm2/s)
                    T Ts = temperature_celsius; // hypothesis T leaf surface = T air
                    T G = street->GetSolarRadiation();
                    T Rsto = ComputeWeselyStomatalResistance(gas_phase_diffusivity_(s), DmH2O, Rsmin_tree, G, Ts);
                    Cut(Rsto);

                    // Compute tree leaves surface resistance
                    T Rs_tree = 1. / ((1. / (Rmes + Rsto)) + (1. / Rcut));
                    Cut(Rs_tree);
                    // Compute dry deposition velocity on tree leaves
                    tree_dry_deposition_velocity = 1. / (Rb_tree + Rs_tree);
                    street->SetTreeDryDepositionVelocity(tree_dry_deposition_velocity,
                                                     DepositionVelocityGlobalIndex(s));
                  }
              }
            else // if the Tree option is activated but the street does not contain any tree
              {
                for (int s = 0; s < Ns_dep; s++)
                  {
                    street->SetTreeDryDepositionVelocity(0.0,
                                                         DepositionVelocityGlobalIndex(s));
                  }
              }
          }
      }
  }

  //! Compute scavening for gaseous species for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::ComputeScavengingCoefficient()
  {
    // Henry constants of scavenged species.
    Array<T, 1> henry_constant_(Ns_scav);
    // Gas phase diffusivities of scavenged species.
    Array<T, 1> gas_phase_diffusivity_(Ns_scav);

    henry_constant_ = 0.0;
    gas_phase_diffusivity_ = 0.0;
    
    // Extracts species data for scavenged species only.
    for (int s = 0; s < Ns_scav; s++)
      {
        henry_constant_(s) = henry_constant[ScavengingName(s)];
        gas_phase_diffusivity_(s) = gas_phase_diffusivity[ScavengingName(s)];
        if (gas_phase_diffusivity_(s) == 0.0)
          throw string("Diffusivity is not given for ") +
            ScavengingName(s);

        if (henry_constant_(s) == 0.0)
          throw string("Henry constant is not given for ") +
            ScavengingName(s);
      }

    Array<T, 1> scavenging_coefficient(Ns_scav);
    scavenging_coefficient = 0.0;
     
    int one = 1;
    Array<T, 1> level(1);
    level = 2.;

    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        T temperature_ = street->GetTemperature();
        T pressure_ = street->GetPressure();
        T rain_ = street->GetRain();
        T cloudheight_ = street->GetHeight();
        T level_ = street->GetHeight();
        int one = 1;

        _compute_scavenging_coefficient(&one, &one, &one, &Ns_scav,
                                        &one, &one, &one,
                                        level.data(),
                                        gas_phase_diffusivity_.data(),
                                        henry_constant_.data(),
                                        &temperature_,
                                        &pressure_,
                                        &rain_,
                                        &cloudheight_,
                                        scavenging_coefficient.data());

        for (int s = 0; s < Ns_scav; s++)
          {
            street->
              SetStreetScavengingCoefficient(scavenging_coefficient(s),
                                             ScavengingGlobalIndex(s));
          }

      }

  }

  //! Sets the scavening global index for a determined specie
  template<class T>
  int StreetNetworkTransport<T>::ScavengingGlobalIndex(int s) const
  {
    return this->GetSpeciesIndex(ScavengingName(s));
  }

  //!Returns the name of the specie in scavening process. 
  template<class T>
  string StreetNetworkTransport<T>::ScavengingName(int s) const
  {
    return species_list_scav.at(s);
  }

  //! Sets the concentrations for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetStreetDryDepositionVelocity()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
          StreetDryDepositionVelocity(s, ist) = street->GetStreetDryDepositionVelocity(s);
        ++ist;
      }
  }

  //! Sets the concentrations for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetWallDryDepositionVelocity()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
          WallDryDepositionVelocity(s, ist) = street->GetWallDryDepositionVelocity(s);
        ++ist;
      }
  }


  //! Sets the concentrations for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetStreetConcentration()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        for (int s = 0; s < this->Ns; ++s)
          StreetConcentration(s, ist) = street->GetStreetConcentration(s);

        ++ist;
      }
  }

  //! Sets the background concentrations for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>
  ::SetStreetBackgroundConcentration(int street_index,
                                     Array<T, 1> background_concentration)
  {
    SetCurrentStreet(street_index);
    Street<T>* street = *current_street;

    for (int s = 0; s < this->Ns; ++s)
        street->SetBackgroundConcentration(background_concentration(s), s);
  }
 
  //! Returns the concentrations for the whole street-network.
  template<class T>
  inline Data<T, 2>& StreetNetworkTransport<T>::GetStreetDryDepositionVelocity()
  {
    return StreetDryDepositionVelocity;
  }

  //! Returns the concentrations for the whole street-network.
  template<class T>
  inline Data<T, 2>& StreetNetworkTransport<T>::GetWallDryDepositionVelocity()
  {
    return WallDryDepositionVelocity;
  }


  //! Returns the concentrations for the whole street-network.
  template<class T>
  inline Data<T, 2>& StreetNetworkTransport<T>::GetStreetConcentration()
  {
    return StreetConcentration;
  }

  //! Returns the number of streets.
  template<class T>
  inline int StreetNetworkTransport<T>::GetNStreet() const
  {
    return total_nstreet;
  }

  //! Returns the number of intersections.
  template<class T>
  inline int StreetNetworkTransport<T>::GetNumberIntersection() const
  {
    return nintersection;
  }

  //! Returns the center coordinate of streets.
  template<class T>
  inline void StreetNetworkTransport<T>::GetStreetCoordinate(int street_index,
                                                             T& longitude,
                                                             T& latitude)
  {
    SetCurrentStreet(street_index);
    Street<T>* street = *current_street;
    longitude = street->GetLongitude();
    latitude = street->GetLatitude();
  }

  //! Returns the center coordinate of intersection.
  template<class T>
  inline void StreetNetworkTransport<T>::
  GetIntersectionCoordinate(int intersection_index,
                            T& longitude,
                            T& latitude)
  {
    SetCurrentIntersection(intersection_index);
    Intersection<T>* intersection = *current_intersection;
    longitude = intersection->GetX();
    latitude = intersection->GetY();
  }

  //! Returns the intersection ID.
  template<class T>
  inline int StreetNetworkTransport<T>::GetIntersectionID(int index)
  {
    SetCurrentIntersection(index);
    Intersection<T>* intersection = *current_intersection;
    return intersection->GetID();
  }

  //! Sets the meteo data for the streets.
  template<class T>
  inline void StreetNetworkTransport<T>::SetStreetMeteo(int street_index,
                                                        T wind_direction,
                                                        T wind_speed, T pblh,
                                                        T ust, T lmo)
  {
    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        if (st == street_index)
          {
            Street<T>* street = *iter;
            street->SetMeteo(wind_direction, wind_speed,
                             pblh, ust, lmo);
          }
        ++st;
      }
  }

  //! Sets the meteo data for the intersections.
  template<class T>
  inline void StreetNetworkTransport<T>::
  SetIntersectionMeteo(int intersection_index,
                       T wind_direction,
                       T wind_speed, T pblh,
                       T ust, T lmo)
  {
    int inter = 0;
    for (typename vector<Intersection<T>* >::iterator iter = IntersectionVector.begin();
         iter != IntersectionVector.end(); iter++)
      {
        if (inter == intersection_index)
          {
            Intersection<T>* intersection = *iter;
            intersection->SetMeteo(wind_direction, wind_speed,
                                   pblh, ust, lmo);
          }
        ++inter;
      }
  }

  //! Compute the concentrations in the street-canyon using the flux balance equation.
  template<class T>
  void StreetNetworkTransport<T>::ComputeStreetConcentration()
  {
    const T delta_concentration_min = 0.01;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        T transfer_velocity = street->GetTransferVelocity(); // m/s
	T temp = transfer_velocity * street->GetWidth() * street->GetLength(); // m3/s
        T outgoing_flux = street->GetOutgoingFlux(); // m3/s
        T street_volume = street->GetHeight() *
          street->GetWidth() * street->GetLength(); // m3
	T street_area = street->GetWidth() * street->GetLength(); // m2
	//! symmetric walls
        T wall_area = 2.0 * street->GetHeight() * street->GetLength(); // m2
        
        bool is_stationary_local = true;
        for (int s = 0; s < this->Ns; ++s)
          {
            T street_conc = street->GetStreetConcentration(s); // ug/m3
            
            T init_street_conc = street->GetInitialStreetConcentration(s); // ug/m3
	    T emission_rate = street->GetEmission(s); //  ug/s
            T inflow_rate = street->GetInflowRate(s); // ug/s
            T conc_bg = street->GetBackgroundConcentration(s); // ug/m3

	    T deposition_rate = 0.0;
	    T scavenging_rate = 0.0;
	    if (this->option_process["with_deposition"])
	      {
	    	T street_dry_deposition_rate = street_area * 
	    	  street->GetStreetDryDepositionVelocity(s); // m3/s
	    	T wall_dry_deposition_rate = wall_area * 
	    	  street->GetWallDryDepositionVelocity(s); // m3/s

                T tree_dry_deposition_rate = 0.0;
                if (this->option_process["with_tree_deposition"])
                  {
                    tree_dry_deposition_rate = street->GetTreeLAI() * street_area *
                      street->GetTreeDryDepositionVelocity(s); // m3/s
                  }
                deposition_rate = street_dry_deposition_rate +
                  wall_dry_deposition_rate + tree_dry_deposition_rate; // m3/s
	      }

	    if (this->option_process["with_scavenging"])
	      {		
		T rain = street->GetRain();
		if (rain > 0.0)
		  scavenging_rate = street_volume * 
		    street->GetStreetScavengingCoefficient(s); // m3/s
	      }

            //! Compute the new concentrations.
            T street_conc_new = 0.0;
	    
            if ((temp + outgoing_flux + deposition_rate + scavenging_rate) == 0.0)
              {
                //! The total outgoing flux should not be zero
                //! because it leads to an infinite increase of pollutant concentrations
                //! in Munich which uses a steady-state approximation.
                //! In this case, the concentrations are assumed to be constant
                //! because the total incoming flux should be equal to zero.
                street_conc_new = (emission_rate * this->Delta_t) / street_volume
                  + init_street_conc;
              }
            else
              street_conc_new = (emission_rate + inflow_rate + temp * conc_bg) /
                (temp + outgoing_flux + deposition_rate + scavenging_rate);

            //! Set the minimum possible concentrations.
            street_conc_new = max(street_conc_new, 0.0);

            //! Check the stationarity
            T delta_concentration = abs(street_conc_new - street_conc);
            if ((street_conc != 0.0) and 
                (delta_concentration > delta_concentration_min))
              is_stationary_local = false;

            //! Set the new concentrations.
            street->SetStreetConcentration(street_conc_new, s);
	    T rain = street->GetRain();
	    T street_dry_deposition_rate = street_area * 
	      street->GetStreetDryDepositionVelocity(s); // m3/s
	    
            T massflux_roof_to_background;
            if (is_stationary_local)
              {
                massflux_roof_to_background = temp * (street_conc_new - conc_bg); // ug/s

                street->SetMassfluxRoofToBackground(massflux_roof_to_background, s);

                T conc_delta = street_conc_new - init_street_conc;
                T street_quantity_delta = conc_delta * street_volume; // ug
                street->SetStreetQuantityDelta(street_quantity_delta, s);
              }

	    if (this->option_process["with_deposition"])
	      street->SetStreetDryDepositionFlux(street_conc_new * 
						 street->GetStreetDryDepositionVelocity(s),
						 s); // ug/m2/s
          }

        //! True if stationarity is obtained for a street.
        street->SetStationary(is_stationary_local);
      }
  }

  //! Compute the concentrations in the street-canyon using the flux balance equation.
  template<class T>
  void StreetNetworkTransport<T>::ComputeStreetConcentrationNoStationary()
  {
    int zero = 0;
    T sub_delta_t_init, sub_delta_t;

    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;

	Array<T, 1> concentration_array(this->Ns);
	Array<T, 1> concentration_array_tmp(this->Ns);
	Array<T, 1> init_concentration_array(this->Ns);
	Array<T, 1> background_concentration_array(this->Ns);
	Array<T, 1> new_concentration_array(this->Ns);
	Array<T, 1> emission_rate_array(this->Ns);
	Array<T, 1> inflow_rate_array(this->Ns);
	Array<T, 1> deposition_rate_array(this->Ns);
	Array<T, 1> street_deposition_rate_array(this->Ns);
	Array<T, 1> scavenging_rate_array(this->Ns);
	Array<T, 1> street_surface_deposited_mass_array(this->Ns);
	Array<T, 1> new_street_surface_deposited_mass_array(this->Ns);
	Array<T, 1> washoff_factor_array(this->Ns);
	
	concentration_array = 0.0;
	concentration_array_tmp = 0.0;
	init_concentration_array = 0.0;
	background_concentration_array = 0.0;
	new_concentration_array = 0.0;
	emission_rate_array = 0.0;
	inflow_rate_array = 0.0;
	deposition_rate_array = 0.0;
	street_deposition_rate_array = 0.0;
	scavenging_rate_array = 0.0;
	street_surface_deposited_mass_array = 0.0;
	new_street_surface_deposited_mass_array = 0.0;
	washoff_factor_array = 0.0;
	
        T transfer_velocity = street->GetTransferVelocity(); // m/s
	T temp = transfer_velocity * street->GetWidth() * street->GetLength(); // m3/s
	T outgoing_flux = street->GetOutgoingFlux(); // m3/s
        T street_volume = street->GetHeight() *
          street->GetWidth() * street->GetLength(); // m3
	T inflow_flux = street->GetIncomingFlux();
        T street_area = street->GetWidth() * street->GetLength(); // m2
	//! symmetric walls
        T wall_area = 2.0 * street->GetHeight() * street->GetLength(); // m2

	for (int s = 0; s < this->Ns; ++s)
          {
	    concentration_array(s) = street->GetStreetConcentration(s);
	    init_concentration_array(s) = street->GetStreetConcentration(s);
	    background_concentration_array(s) = street->GetBackgroundConcentration(s);
	    emission_rate_array(s) = street->GetEmission(s);
	    inflow_rate_array(s) = street->GetInflowRate(s);
	    if (this->option_process["with_deposition"])
	      {
                // Deposition flux on the street ground
		street_deposition_rate_array(s) = street_area * 
	    	  street->GetStreetDryDepositionVelocity(s); // m3/s

	    	T wall_dry_deposition_rate = wall_area * 
	    	  street->GetWallDryDepositionVelocity(s); // m3/s

                T tree_dry_deposition_rate = 0.0;
                if (this->option_process["with_tree_deposition"])
                  {
                    tree_dry_deposition_rate = street->GetTreeLAI() * street_area *
                      street->GetTreeDryDepositionVelocity(s); // m3/s
                  }

                // Total deposition flux on the ground, walls and tree leaves of the street
	    	deposition_rate_array(s) = street_deposition_rate_array(s) +
                  wall_dry_deposition_rate + tree_dry_deposition_rate; // m3/s
	      }
	    if (this->option_process["with_scavenging"])
	      {
		T rain = street->GetRain();
		if (rain > 0.0)
		  scavenging_rate_array(s) = street_volume * 
		    street->GetStreetScavengingCoefficient(s); // m3/s
	      }
	  }

	InitStep(sub_delta_t_init,
		 sub_delta_t_min,
		 transfer_velocity,
		 temp,
		 outgoing_flux,
		 street_volume,
		 concentration_array,
		 background_concentration_array,
		 emission_rate_array,
		 inflow_rate_array,
		 deposition_rate_array,
		 street_deposition_rate_array,
		 scavenging_rate_array,
		 zero, //resuspension_factor
		 washoff_factor_array,
		 street_surface_deposited_mass_array, //zero for gas species
		 this->Ns,
		 zero);

	Date current_date_tmp = this->current_date;
	Date next_date = this->current_date;
	Date next_date_tmp = this->current_date;
	next_date.AddSeconds(this->Delta_t);
	next_date_tmp.AddSeconds(sub_delta_t_init);
	
	while (current_date_tmp < next_date)
	  {
	    //! Get street concentrations.
	    for (int s = 0; s < this->Ns; ++s)
	      concentration_array(s) = street->GetStreetConcentration(s);

            //! Use the ETR method to calculates new street concentrations.
	    if (option_method == "ETR")
	      {
		ETRConcentration(transfer_velocity,
				 temp,
				 outgoing_flux,
				 street_volume,
				 concentration_array,
				 concentration_array_tmp,
				 background_concentration_array,
				 emission_rate_array,
				 inflow_rate_array,
				 deposition_rate_array,
				 street_deposition_rate_array,
				 scavenging_rate_array,
				 zero, //resuspension_factor
				 washoff_factor_array,
				 street_surface_deposited_mass_array, //zero for gas species
				 new_concentration_array,
				 sub_delta_t_init,
				 this->Ns,
				 zero,
				 new_street_surface_deposited_mass_array);
	      }
	    else if(option_method == "Rosenbrock")
	      {
		RosenbrockConcentration(transfer_velocity,
					temp,
					outgoing_flux,
					street_volume,
					concentration_array,
					concentration_array_tmp,
					background_concentration_array,
					emission_rate_array,
					inflow_rate_array,
					deposition_rate_array,
					new_concentration_array,
					sub_delta_t_init,
					inflow_flux,
					street_deposition_rate_array,
					scavenging_rate_array,
					zero, //resuspension_factor
					washoff_factor_array,
					street_surface_deposited_mass_array, //zero for gas species
					new_street_surface_deposited_mass_array,
					this->Ns,
					zero);
	      }
	    else
	      throw string("Error: numerical method not chosen.");
	    
	    
	    //! Set the new concentrations.
	    for (int s = 0; s < this->Ns; ++s)
	      street->SetStreetConcentration(new_concentration_array(s), s);

	    //! Calculates the new sub_delta_t for the next iteration
	    T sub_delta_t_max = next_date.GetSecondsFrom(next_date_tmp); 
	    sub_delta_t_max = max(sub_delta_t_max, 0.0);
	    AdaptTimeStep(new_concentration_array,
			  concentration_array_tmp,
			  sub_delta_t_init,
			  sub_delta_t_min,
			  sub_delta_t_max,
			  sub_delta_t);

	    //! Actualises current_time_tmp
	    current_date_tmp.AddSeconds(sub_delta_t_init);
	    next_date_tmp.AddSeconds(sub_delta_t);
	    
	    //! Set the new sub_delta_t
	    sub_delta_t_init = sub_delta_t;
	  }
	
	for (int s = 0; s < this->Ns; ++s)
	  {
            T massflux_roof_to_background;
	    massflux_roof_to_background = temp * (new_concentration_array(s) - background_concentration_array(s)); // ug/s
	    
	    street->SetMassfluxRoofToBackground(massflux_roof_to_background, s);
   
	    T conc_delta = new_concentration_array(s) - init_concentration_array(s);
	    T street_quantity_delta = conc_delta * street_volume; // ug
	    street->SetStreetQuantityDelta(street_quantity_delta, s);
	  }
      }
  }
  //! Compute the concentrations in the street-canyon using the flux balance equation.
  template<class T>
  void StreetNetworkTransport<T>::SetInitialStreetConcentration()
  {
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
          {
            T street_conc = street->GetStreetConcentration(s); // ug/m3
            street->SetInitialStreetConcentration(street_conc, s);
          }
      }
  }


  //! Compute the entering/outgoing flow rate for the intersection.
  template<class T>
  void StreetNetworkTransport<T>::ComputeIntersectionFlux(Array<T, 2>& extended_matrix,
                                                          T wind_dir_inter)
  {
    int i, j;

    int nstreet_inter = StreetVectorInter.size(); 

    // Compute wind speed and flux in the streets
    Array<double, 1> flux(nstreet_inter); // Volume flux in m3/s
    Array<double, 1> temp(nstreet_inter);

    // Define wind directions in the streets: upwind or downwind
    int nstreet_in = 0; // Number of incoming streets
    int nstreet_out = 0; // Number of outgoing streets
    Array<int, 1> ind_wind(nstreet_inter);

    for (int i = 0; i < nstreet_inter; i++)
      {
        Street<T>* street = StreetVectorInter[i];
        flux(i) = street->GetStreetWindSpeed() * street->GetWidth() * street->GetHeight();
        temp(i) = street->GetStreetAngleIntersection() * 180. / pi;

        T dangle = abs(street->GetStreetAngleIntersection() - wind_dir_inter);
        
        if (dangle > (pi / 2.0) and dangle < (pi * 3.0 / 2.0))
          {
            ind_wind(i) = 0; // incoming flow
            nstreet_in += 1;
          }
        else
          {
            ind_wind(i) = 1; // outgoing flow
            nstreet_out += 1;
          }
      }

    //! Additional arrays
    Array<T, 2> alpha(nstreet_in, nstreet_out); // Flux ratio 
    alpha = 0.0;
    Array<T, 1> P_in(nstreet_in), P_out(nstreet_out); // Flux 

    //! Set the reference incoming street: first upwind(incoming) street in the counterclockwise.
    Array<int, 1> ind_street_in(nstreet_in);

    //! Get the incoming street numbers
    j = 0;
    for (int i = 0; i < nstreet_inter; i++)
      if (ind_wind(i) == 0)
        {
          ind_street_in(j) = i;
          ++j;
        }

    // Sort in the order of decreasing angle
    for (int i = 0; i < nstreet_in; i++)
      for (int j = i; j < nstreet_in; j++)
      {
        int a1 = ind_street_in(i);
        int a2 = ind_street_in(j);
        if (StreetVectorInter[a1]->GetStreetAngleIntersection() < 
            StreetVectorInter[a2]->GetStreetAngleIntersection())
          {
            ind_street_in(i) = a2;
            ind_street_in(j) = a1;
          }
      }

    //! Sort in the counterclockwise order
    bool sorted = false;
    while ((nstreet_in > 1) and (!sorted))
      {
        if (((StreetVectorInter[ind_street_in(0)]->GetStreetAngleIntersection()) - (StreetVectorInter[ind_street_in(1)]->GetStreetAngleIntersection())) <= pi)
        sorted = true;
      else
        {
          int temp = ind_street_in(0);
          for (int i = 0; i < (nstreet_in - 1); i++)
            ind_street_in(i) = ind_street_in(i + 1);
          ind_street_in(nstreet_in - 1) = temp;
        }
      }

    //! Set the reference outgoing street: first downwind(outgoing) street in the clockwise.
    Array<int, 1> ind_street_out(nstreet_out);
    //! Get the outgoing street numbers
    j = 0;
    for (int i = 0; i < nstreet_inter; i++)
      if (ind_wind(i) == 1)
        {
          ind_street_out(j) = i;
          ++j;
        }

    //! Sort in increasing order
    for (int i = 0; i < nstreet_out; i++)
      for (int j = i; j < nstreet_out; j++)
        {
          int a1 = ind_street_out(i);
          int a2 = ind_street_out(j);
          if (StreetVectorInter[a1]->GetStreetAngleIntersection() > StreetVectorInter[a2]->GetStreetAngleIntersection())
            {
              ind_street_out(i) = a2;
              ind_street_out(j) = a1;
            }
        }

    //! Sort in clockwise order
    sorted = false;
    while ((nstreet_out > 1) and (!sorted))
      {
        if (((StreetVectorInter[ind_street_out(nstreet_out - 1)]->
              GetStreetAngleIntersection()) - 
             (StreetVectorInter[ind_street_out(nstreet_out - 2)]->
              GetStreetAngleIntersection())) <= pi)
          sorted = true;
        else
          {
            int temp = ind_street_out(nstreet_out - 1);
            for (int i = 0; i < (nstreet_out - 1); i++)
              ind_street_out(nstreet_out - 1 - i) = ind_street_out(nstreet_out - 2 - i);
            ind_street_out(0) = temp;
          }
      }

    //! Set incoming/outgoing flux
    for (int i = 0; i < nstreet_in; i++)
      P_in(i) = flux(ind_street_in(i));
    for (int i = 0; i < nstreet_out; i++)
      P_out(i) = flux(ind_street_out(i));

    //! Compute flux exchange with the atmosphere
    T P0 = 0.0; // Flux of exchange with the atmosphere
    //! P0 > 0: outgoing flux from the intersection
    //! P0 < 0: incoming flux into the intersection
    T sum_P_in = 0;
    for (int i = 0; i < nstreet_in; i++)
      sum_P_in += P_in(i);
    T sum_P_out = 0;
    for (int i = 0; i < nstreet_out; i++)
      sum_P_out += P_out(i);
    P0 = sum_P_in - sum_P_out;

    //! Distribution of P0
    Array<T, 2> flux_matrix(nstreet_in + 1, nstreet_out + 1); 
    flux_matrix = 0.0;
    T alpha0 = 0;
    if (P0 > 0) // outgoing flux
      {
        alpha0 = P0 / sum_P_in;
        for (int i = 0; i < nstreet_in; i++)
          {
            flux_matrix(i + 1, 0) = alpha0 * P_in(i); 
            P_in(i) = (1 - alpha0) * P_in(i);
          }
      }
    else if (P0 < 0) // incoming flux
      {
        alpha0 = abs(P0) / sum_P_out;
        for (int i = 0; i < nstreet_out; i++)
          {
            flux_matrix(0, i + 1) = alpha0 * P_out(i); 
            P_out(i) = (1 - alpha0) * P_out(i);
          }
      }

    sum_P_in = 0.0;
    sum_P_out = 0.0;
    for (i = 0; i < nstreet_in; i++)
      sum_P_in += P_in(i);
    for (i = 0; i < nstreet_out; i++)
      sum_P_out += P_out(i);

    //! Compute flux among the streets.
    ComputeAlpha(nstreet_in, nstreet_out, P_in, P_out, alpha, flux_matrix);

    CreateExtendedMatrix(ind_street_in, ind_street_out, flux_matrix, 
                         extended_matrix);
  }


  //! Compute fluxes at the intersection using the gaussian distribution
  //! in the horizontal wind direction.
  /*!
    \param gaussian_factor 
    \param intersection intersection object
  */
  template<class T>
  void StreetNetworkTransport<T>::ComputeGaussianFluxMatrix(T gaussian_factor, 
                                                            Intersection<T>* intersection)
  {
    int nstreet_inter = intersection->GetNStreet();
    Array<T, 2> flux_matrix(nstreet_inter + 1, nstreet_inter + 1);
    Array<T, 2> gaussian_matrix(nstreet_inter + 1, nstreet_inter + 1);
    flux_matrix = intersection->GetFluxMatrix();
    gaussian_matrix = intersection->GetGaussianMatrix();
    gaussian_matrix += flux_matrix * gaussian_factor;
    intersection->SetGaussianMatrix(gaussian_matrix);
  }

  //! Compute the street angle which is defined as the angle 
  //! bewteen the intersections which the street connects 
  //! (begin intersection to end intersection).
  /*!
    \param ind_street_in
    \param ind_street_out
    \param flux_matrix
    \param extended_matrix
  */
  template<class T>
  void StreetNetworkTransport<T>::CreateExtendedMatrix(Array<int, 1> ind_street_in,
                                                       Array<int, 1> ind_street_out,
                                                       Array<T, 2> flux_matrix,
                                                       Array<T, 2>& extended_matrix)
  {
    extended_matrix = 0.0;
    for (int i = 0; i < flux_matrix.shape()[0]; i++)
      for (int j = 0; j < flux_matrix.shape()[1]; j++)
        {
          if (i == 0 and j == 0)
            extended_matrix(i, j) = 0.0;
          else if (i == 0 and j != 0)
            extended_matrix(i, ind_street_out(j - 1) + 1) = flux_matrix(i, j);
          else if (j == 0 and i != 0)
            extended_matrix(ind_street_in(i - 1) + 1, j) = flux_matrix(i, j);
          else
            extended_matrix(ind_street_in(i - 1) + 1, ind_street_out(j - 1) + 1) 
              = flux_matrix(i, j);
        }
  }


  //! Compute the street angle which is defined as the angle bewteen the intersections which the street connects (begin intersection to end intersection).
  template<class T>
  void StreetNetworkTransport<T>::ComputeStreetAngle()
  {
    T x1, x2, y1, y2;
    x1 = 0.0;
    x2 = 0.0;
    y1 = 0.0;
    y2 = 0.0;
    bool begin_matched = false;
    bool end_matched = false;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        for (typename vector<Intersection<T>* >::iterator iter2 = IntersectionVector.begin(); iter2 != IntersectionVector.end(); iter2++)
          {
            Intersection<T>* intersection = *iter2;
            if (street->GetBeginIntersectionID() == intersection->GetID())
              {
                x1 = intersection->GetX();
                y1 = intersection->GetY();
                begin_matched = true;
              }
            if (street->GetEndIntersectionID() == intersection->GetID())
              {
                x2 = intersection->GetX();
                y2 = intersection->GetY();
                end_matched = true;
              }
          }
        if ((!begin_matched) or (!end_matched))
          throw string("Wrong intersection index is given.");      

        //! Set the coordinate
        street->SetCoordinate((x1 + x2) * 0.5, (y1 + y2) * 0.5);

        T dx = x2 - x1;
        T x_distance = abs(earth_radius * dx * pi / 180.);
        T dy = y2 - y1;
        T y_distance = abs(earth_radius * (sin(y2 * pi / 180.) - 
                                                sin(y1 * pi / 180.)));
        T dl = sqrt(pow(x_distance, 2) + pow(y_distance, 2));
        if (dl == 0.0)
          throw string("Distance between the intersections is zero.");

        //! gamma is arc cosine of input value, x.
        //! 0 < x < 1 and 0 < gamma < PI/2
        //! When x is 1, gamma is 0 radian.
        //! When x is 0, gamma is PI/2 radian.
        T gamma = acos(abs(x_distance / dl));

        //! ********** Street angle **********
        //! Begin intersection: the intersection 1.
        //! End intersection: the intersection 2.
        //! When the intersection 2 is located north of the intersectoin 1, 
        //! the street angle is defined as zero. And it increases clockwise from north.
        //! **********************************

        //! Case 1: the intersection 2 is located northeast of the intersection 1
        //! (0 =< street_angle =< PI/2 in radian)
        //! When x is 1 (x_distance is equal to dl), gamma is 0 radian.
        //! ---> It represents east and the street angle is PI/2 clockwise from north. 
        //! When x is 0 (x_distance is 0), gamma is PI/2 radian.
        //! ---> It represents north and the street angle is 0 from north. 
        if ((dx >= 0.0) and (dy >= 0.0))
          street->SetStreetAngle(pi / 2.0 - gamma);
        //! Case 2: the intersection 2 is located southeast of the intersection 1
        //! (PI/2 < street_angle =< PI in radian)
        //! When x is 1 (x_distance is equal to dl), gamma is 0 radian.
        //! ---> It represents east and the street angle is PI/2 clockwise from north. 
        //! When x is 0 (x_distance is 0), gamma is PI/2 radian.
        //! ---> It represents south and the street angle is PI clockwise from north. 
        else if ((dx >= 0.0) and (dy < 0.0))
          street->SetStreetAngle(pi / 2.0 + gamma);
        //! Case 3: the intersection 2 is located southwest of the intersection 1
        //! (PI < street_angle =< PI*3/2 in radian)
        //! When x is 1 (x_distance is equal to dl), gamma is 0 radian.
        //! ---> It represents west and the street angle is PI*3/2 clockwise from north. 
        //! When x is 0 (x_distance is 0), gamma is PI/2 radian.
        //! ---> It represents south and the street angle is PI clockwise from north. 
        else if ((dx < 0.0) and (dy <= 0.0))
          street->SetStreetAngle(pi * 1.5 - gamma);
        //! Case 4: the intersection 2 is located northwest of the intersection 1
        //! (PI*3/2 < street_angle =< PI*2 in radian)
        //! When x is 1 (x_distance is equal to dl), gamma is 0 radian.
        //! ---> It represents west and the street angle is PI*3/2 clockwise from north. 
        //! When x is 0 (x_distance is 0), gamma is PI/2 radian.
        //! ---> It represents north and the street angle is PI*2 clockwise from north. 
        else if ((dx < 0.0) and (dy > 0.0))
          street->SetStreetAngle(pi * 1.5 + gamma);
      }
  }

  //! Compute the factor for the Gaussian distribution.
  template<class T>
  T StreetNetworkTransport<T>::ComputeGaussian(double theta, double sigma_theta, 
                                      double theta0)
  {
    double gaussian = 1.0 / (sigma_theta * sqrt(2 * pi)) * exp(-0.5 * pow((theta - theta0) / sigma_theta, 2));
    return gaussian;
  }

  //! Compute the sigma_v values based on the grid-averaged friction velocity.
  //! \param ustar grid-averaged friction velocity (m/s)
  //! Based on Hunt et al. (1988), CERC (2001), and Souhlac et al. (2011)
  template<class T>
  T StreetNetworkTransport<T>::ComputeSigmaV(T lmo, T pblh, T ust)
  {
    int nz = 10;
    T sigma_v_ = 0.0;
    Array<T, 1> z(nz);
    for (int k = 0; k < nz; ++k)
      {
        z(k) = pblh / (nz - 1) * k;
        if (lmo < 0.0)
          {
            T wst = ust * pow((pblh / (karman * abs(lmo))), (1.0 / 3.0));
            sigma_v_ += sqrt(0.3 * pow(wst, 2.0) + 
                             pow((2.0 * ust * (1.0 - 0.8 * z(k) / pblh)), 2.0));
          }
        else if (lmo >= 0.0 and lmo < pblh)
          sigma_v_ += 2.0 * ust * pow((1.0 - 0.5 * z(k) / pblh), 3.0 / 4.0);
        else
          sigma_v_ += 2.0 * ust * (1.0 - 0.8 * z(k) / pblh);
      }
    sigma_v_ /= nz;

    return sigma_v_; 
  }

  //! Compute the standard deviation of the vertical wind speed at a roof.
  //! Based on Hunt et al. (1988), CERC (2001), and Soulhac et al. (2011).
  template<class T>
  void StreetNetworkTransport<T>::ComputeSigmaW()
  {
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        T sigma_w_ = 0.0;
        T ust = street->GetStreetUstar();
        T pblh = street->GetPBLH();
        T lmo = street->GetLMO();
        T height = street->GetHeight();

        //! Unstable
        if (lmo < 0.0)
          {
            T wst = ust * pow((pblh / (karman * abs(lmo))), (1.0 / 3.0));
            T sigma_wc = sqrt(0.4) * wst * 2.1 * pow((height / pblh), (1.0 / 3.0))
              * (1.0 - 0.8 * height / pblh);
            T sigma_wn = 1.3 * ust * (1.0 - 0.8 * height / pblh);
            sigma_w_ = sqrt(pow(sigma_wc, 2.0) + pow(sigma_wn, 2.0));
          }
        //! Stable
        else if (lmo >= 0.0 and lmo < pblh)
          sigma_w_ = 1.3 * ust * pow((1.0 - 0.5 * height / pblh), 3.0 / 4.0);
        //! Neutral
        else
          sigma_w_ = 1.3 * ust * (1.0 - 0.8 * height / pblh);

        if (sigma_w_ != sigma_w_)
          throw string("Error: Nan in sigma_w");
        street->SetSigmaW(sigma_w_);
      }
  }

  //! Compute the transfer velocity at roof level.
  template<class T>
  void StreetNetworkTransport<T>::ComputeTransferVelocity()
  {
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        T H = street->GetHeight();
        T W = street->GetWidth();
        T sigma_w = street->GetSigmaW();
        T min_velocity = 0.001;
        T velocity = 0.0;
        if (option_transfer == "Sirane")
          velocity = sigma_w / (sqrt(2.0) * pi);
        else if (option_transfer == "Schulte")
          {
            T lm = ComputeSchulteLm(H, W);
            velocity = sigma_w * lm;
          }
         else if (option_transfer == "Wang")
          {
            T sH;
            if (this->option_process["with_tree_aerodynamic"])
              sH = ComputeWangsH(H, W, street->GetTreeHeight(), street->GetTreeLAI(), Cdt);
            else
              sH = ComputeWangsH(H, W, 0.0, 0.0, 0.0);
            velocity = sigma_w * karman * sH;
          }
        street->SetTransferVelocity(max(velocity, min_velocity));
      }
  }


  //! Compute ustar and uH with the method of Macdonald et al (1998)
  template<class T>
  void StreetNetworkTransport<T>::ComputeMacdonaldProfile()
  {

    T LambdaP = this->building_density; // plan area density of obstacles
    T building_width = 0.5 * LambdaP / (1.0 - LambdaP) * Mean_width;
    T Af = Mean_length * Mean_height; // frontal area of obstacles
    T At = Mean_length * (Mean_width + building_width); // lot area of obstacles
    T LambdaF = Af / At; // frontal area density of obstacles

    T Cd_building = 1.2; // Drag coefficient for buildings
    T A = 4.43; // Correction coefficient A = 4.43 for staggered arrays and A = 3.59 for square array (Macdonald et al, 1998)
    T b = 1.0; // Correction coefficient b = 1.0 for staggered arrays and b = 0.55 for square array (Macdonald et al, 1998)

    d_city = ComputeMacdonaldD(Mean_height, LambdaP, A);
    z0_city = ComputeMacdonaldZ0(Mean_height, LambdaF, LambdaP, d_city, b, Cd_building);

    T ustar_macd;

    if (option_uH_macdonald == "Uref")
      {
        for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
          {
            Street<T>* street = *iter;
            // Compute average Macdonald ustar on street network
            T u_zref = street->GetWindSpeed(); // wind speed at the reference altitude (m/s)
            ustar_macd = ComputeMacdonaldUstar(u_zref, zref, d_city, z0_city);
            street->SetStreetUstar(ustar_macd);

            // Compute Macdonald uH for each street of the street network
            T H = street->GetHeight();
            if (zref < H)
              throw string("Error: reference altitude is inferior to the building height in the street ") + to_str(street->GetStreetID())+ ".";
          }
    
        // Friction velocity using Macdonald method for intersections in the street network
        for (typename vector<Intersection<T>* >::iterator iter = IntersectionVector.begin();
             iter != IntersectionVector.end(); iter++)
          {
            Intersection<T>* intersection = *iter;
            // wind speed at the reference altitude (m/s)       
            T u_zref = intersection->GetWindSpeed();
            ustar_macd = ComputeMacdonaldUstar(u_zref, zref, d_city, z0_city);
            intersection->SetIntersectionUST(ustar_macd);
          }
      }
  }

  
  //! Compute the wind speed in the street-canyon.
  //! Option "Sirane" based on Soulhac et al. (2008)
  //! Option "Exponential" based on Lemonsu et al. (2004), and Cherin et al. (2015)
  //! Option "Wang" based on Wang (2012; 2014)
  template<class T>
  void StreetNetworkTransport<T>::ComputeUstreet()
  {
    T Ustreet_y;
    T Ustreet = 0.0;
    T Uh = 0.0;
    T Um = 0.0;

    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        T H = street->GetHeight();
        T W = street->GetWidth();       
        T ang = street->GetStreetAngle();
        T wind_direction = street->GetWindDirection();
        T phi = abs(wind_direction - ang);
        T ustar = street->GetStreetUstar();
        int nz = int(10 * H);
        if (option_uH == "Sirane")
          {
            T delta = min(H, W / 2);
            T C = ComputeSiraneC(z0s, delta);
            Um = ComputeSiraneUm(C, ustar, z0s);
            T Uh_Um_factor = ComputeUmtoUHFactor(C, z0s);
            if (option_ustreet == "Sirane")
              Ustreet = ComputeSiraneUstreet(H, W, C, z0s, delta, Um);
            else if (option_ustreet == "Wang")
              {
                T sH;
                Uh = Um * Uh_Um_factor;
                if (this->option_process["with_tree_aerodynamic"])
                  {
                    sH = ComputeWangsH(H, W, street->GetTreeHeight(), street->GetTreeLAI(), Cdt);
                    Ustreet = ComputeWangUstreet(H, W, phi, z0s, sH, nz, Uh, street->GetTreeLAI(), Cdt);
                  }
                else
                  {
                    sH = ComputeWangsH(H, W, 0.0, 0.0, 0.0);
                    Ustreet = ComputeWangUstreet(H, W, phi, z0s, sH, nz, Uh, 0.0, 0.0);
                  }
              }
            else if (option_ustreet == "Exponential")
              {
                Uh = Um * Uh_Um_factor;
                Ustreet = ComputeExpUstreet(H, W, z0s, Uh);
              }
            else
              throw("Wrong option given. Choose Sirane, Wang or Exponential.");
          }
        else if (option_uH == "Macdonald")
          {
            ComputeMacdonaldProfile();
            if  (H < (d_city + z0_city))
              Uh = 0.0;
            else
              Uh = ComputeMacdonaldUH(H, ustar, d_city, z0_city);
            if (option_ustreet == "Sirane")
              {
                T delta = min(H, W / 2);
                T C = ComputeSiraneC(z0s, delta);
                T Uh_Um_factor = ComputeUmtoUHFactor(C, z0s);
                Um = Uh / Uh_Um_factor;
                Ustreet = ComputeSiraneUstreet(H, W, C, z0s, delta, Um);
              }
            else if (option_ustreet == "Wang")
              {
                T sH;
                if (this->option_process["with_tree_aerodynamic"])
                  {
                    sH = ComputeWangsH(H, W, street->GetTreeHeight(), street->GetTreeLAI(), Cdt);
                    Ustreet = ComputeWangUstreet(H, W, phi, z0s, sH, nz, Uh, street->GetTreeLAI(), Cdt);
                  }
                else
                  {
                    sH = ComputeWangsH(H, W, 0.0, 0.0, 0.0);
                    Ustreet = ComputeWangUstreet(H, W, phi, z0s, sH, nz, Uh, 0.0, 0.0);
                  }
              }
            else if (option_ustreet == "Exponential")
              Ustreet = ComputeExpUstreet(H, W, z0s, Uh);
            else
              throw("Wrong option given. Choose Sirane, Wang or Exponential.");
          }
        else
          throw("Wrong option given. Choose Sirane or Macdonald.");

        Ustreet_y = Ustreet * abs(cos(phi));
        street->SetStreetWindSpeed(max(Ustreet_y, ustreet_min));
      }
  }

  //! Compute the horizontal fluctuation of the wind direction.
  template<class T>
  void StreetNetworkTransport<T>::ComputeWindDirectionFluctuation()
  {
    //! Maximum sigma_theta: 10° based on Ben Salem et al. (2015) 
    //! Maximum fluctuation of +/-20° (2 * sigma_theta)
    //! in radian
    const double max_sigma_theta(pi / 18.); 

    for (typename vector<Intersection<T>* >::iterator iter = IntersectionVector.begin(); iter != IntersectionVector.end(); iter++)
      {
        //! Get the intersection class object.
        Intersection<T>* intersection = *iter;
        T theta0 = intersection->GetWindDirection();        

        if (this->option_process["with_horizontal_fluctuation"])
          {
            T sigma_v = ComputeSigmaV(intersection->GetLMO(), intersection->GetPBLH(),
                                      intersection->GetUST());
            //! Blackadar (1997) and Soulhac et al. (2009)
            //! in radian
            T sigma_theta = min(sigma_v / intersection->GetWindSpeed(), max_sigma_theta);
            intersection->SetSigmaTheta(sigma_theta);
            //! Maximum ntheta is 10 when sigma_theta is equal to the maximum value.
            int ntheta = int(sigma_theta / pi * 180.);
            if (ntheta > 1)
              {
                Array<T, 1> theta(ntheta), gaussian(ntheta);
                T step = 4.0 * sigma_theta / (ntheta - 1);
                T begin = (theta0 - 2.0 * sigma_theta);
                for (int i = 0; i < ntheta; i++)
                  {
                    theta(i) =  begin + i * step;
                    gaussian(i) = ComputeGaussian(theta(i), sigma_theta, theta0);
                    if (theta(i) < 0.0)
                      theta(i) += (2.0 * pi);
                    
                    //! Compute the flux at the intersection.
                    ComputeIntersection(theta(i), intersection);

                    ComputeGaussianFluxMatrix(step * gaussian(i), intersection);

                  }
              }
            else
              ComputeIntersection(theta0, intersection);
          }
        else
          ComputeIntersection(theta0, intersection);

      }
  }

  //! Compute the flux distribution rate in the intersection.
  //! This is not indispensable for the model.
  template<class T>
  void StreetNetworkTransport<T>::ComputeAlpha(int nstreet_in, int nstreet_out, 
                               Array<double, 1> P_in, Array<double, 1> P_out,
                               Array<double, 2>& alpha, Array<double, 2>& flux_matrix)
  {
    // General parameters and options 
    int i, j;
    double sum_alpha = 0.0;
    double temp_in;
    Array<double, 1> temp_out(nstreet_out);
    temp_out = P_out;
    for (i = 0; i < nstreet_in; i++)
      {
        temp_in = P_in(i);
        for (j = 0; j < nstreet_out; j++)
          {
            if (temp_in > temp_out(j))
              flux_matrix(i + 1, j + 1) = temp_out(j);
            else
              flux_matrix(i + 1, j + 1) = temp_in;
            temp_in -= flux_matrix(i + 1, j + 1);
            temp_out(j) -= flux_matrix(i + 1, j + 1);
            if (P_in(i) != 0.0)
              alpha(i, j) = flux_matrix(i + 1, j + 1) / P_in(i);
          }
      }
  }


  //! Compute the inflow rate with the extended matrix.
  template<class T>
  void StreetNetworkTransport<T>::ComputeInflowRateExtended()
  {
    for (typename vector<Intersection<T>* >::iterator iter = IntersectionVector.begin();
         iter != IntersectionVector.end(); iter++)
      {
        //! Get the intersection class object.
        Intersection<T>* intersection = *iter;

        //! Get the number of streets that are connected to the intersection 
        //! and their list.
        int nstreet_inter = intersection->GetNStreet(); 
        Array<int, 1> street_list_inter(nstreet_inter);
        street_list_inter = intersection->GetStreetList();

        //! Get the matrix for the flux distribution at the intersection.
        //! ====== Description of the matrix =====
        //! --- First line : Flux from the atmosphere to each street.
        //! --- First column: Flux from each street to the atmosphere.
        //! --- (n, m): Flux from the n_th street to the m_th street.
        //! --- Diagonal flux are zero.
        Array<T, 2> extended_matrix(nstreet_inter + 1, nstreet_inter +1);
        if (this->option_process["with_horizontal_fluctuation"])
          extended_matrix = intersection->GetGaussianMatrix();
        else
          extended_matrix = intersection->GetFluxMatrix();


        for (int j = 0; j < nstreet_inter; ++j)
          {
            int end_inter_id = 0;
            int begin_inter_id = 0;
            for (typename vector<Street<T>* >::iterator iter2 = StreetVector.begin();
                 iter2 != StreetVector.end(); iter2++)
              {
                Street<T>* street = *iter2;
                if (street_list_inter(j) == street->GetStreetID())
                  {
                    StreetVectorInter.push_back(street);
                    break;
                  }
              }
          }

        //! Compute outgoing volume flux (m3/s)
        for (int i = 0; i < nstreet_inter + 1; i++)
          for (int j = 0; j < nstreet_inter + 1; j++)
            if (i != j and i != 0)
              {
                T old_outgoing_flux = StreetVectorInter[i - 1]->GetOutgoingFlux();
                T outgoing_flux = extended_matrix(i, j) + old_outgoing_flux;
                StreetVectorInter[i - 1]->SetOutgoingFlux(outgoing_flux);
              }

	//! Compute incoming volume flux from other streets (m3/s)
	for (int i = 0; i < nstreet_inter + 1; i++)
	  for (int j = 0; j < nstreet_inter + 1; j++)
	    if (i != j and i != 0 and j != 0)
	      {
		T old_incoming_flux = StreetVectorInter[j - 1]->GetIncomingFlux();
		T incoming_flux = extended_matrix(i, j) + old_incoming_flux;
		
		StreetVectorInter[j - 1]->SetIncomingFlux(incoming_flux);
	      }

        for (int s = 0; s < this->Ns; ++s)
          {
            for (int i = 0; i < nstreet_inter + 1; i++)
              for (int j = 0; j < nstreet_inter + 1; j++)
                {
                  if (i != j and i == 0 and j != 0)
                    {
                      T old_massflux = StreetVectorInter[j - 1]->GetMassfluxFromBackground(s);
                      T old_inflow_rate = StreetVectorInter[j - 1]->GetInflowRate(s);

                      //! Compute the inflow rate from the atmosphere to the street
                      T new_massflux  = StreetVectorInter[j - 1]->GetBackgroundConcentration(s) * extended_matrix(i, j);
                      T massflux_from_background = new_massflux + old_massflux;
                      T inflow_rate = new_massflux + old_inflow_rate;
                      StreetVectorInter[j - 1]->SetMassfluxFromBackground(massflux_from_background, s);
                      StreetVectorInter[j - 1]->SetInflowRate(inflow_rate, s);
                    }
                  else if (i != j and i != 0 and j == 0)
                    {
                      T old_massflux = StreetVectorInter[i - 1]->GetMassfluxToBackground(s);
                      //! Compute the outflow rate to the atmosphere from the street
                      T new_massflux = StreetVectorInter[i - 1]->GetStreetConcentration(s) * extended_matrix(i, j);
                      T massflux_to_background = old_massflux + new_massflux;
                      StreetVectorInter[i - 1]->SetMassfluxToBackground(massflux_to_background, s);
                    }
                  else if (i != j and i != 0 and j != 0)
                    {
                      T new_massflux = StreetVectorInter[i - 1]->GetStreetConcentration(s) * extended_matrix(i, j);
                      T old_inflow_rate = StreetVectorInter[j - 1]->GetInflowRate(s);
                      T inflow_rate = new_massflux + old_inflow_rate;
                      StreetVectorInter[j - 1]->SetInflowRate(inflow_rate, s);
                    }
                }
          }
        
        //! Clear the street vector for the intersection. 
        StreetVectorInter.clear();
      }
  }

  //! Write the results in the output files in the text format.
  template<class T>
  void StreetNetworkTransport<T>::OutputSaver()
  {
    if (text_file)
      {
        Date previous_date = this->current_date;
        previous_date.AddSeconds(this->Delta_t * -1.0);

        ofstream output_file, output_file_date;
        string hour, minute;
        ostringstream sh, sm;
        sh << setw(2) << setfill('0') << previous_date.GetHour();
        hour = sh.str();
        sm << setw(2) << setfill('0') << previous_date.GetMinutes();
        minute = sm.str();
        string output_date = to_str(previous_date.GetDate()) + "_" +  
          hour + "-" + minute;

        string file_name_date = output_dir + "/" + output_date;
        output_file_date.open(file_name_date.c_str(), ios_base::trunc);
        for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); 
             iter != StreetVector.end(); iter++)
          {
            Street<T>* street = *iter;
            int id = street->GetStreetID();
            output_file_date << to_str(id) << "\t";
            for (int s = 0; s < this->Ns; ++s)
              {
                // Leighton
                if (this->Ns == 4)
                  {
                    // NO2: s = 2 
                    // NO: s = 3
                    // O3: s = 1
                    if (s == 1 or s == 2 or s == 3) 
                      output_file_date << street->GetStreetConcentration(s) << "\t" ;
                  }
                // CB05
                else
                  {
                    // NO2: s = 51 
                    // NO: s = 47
                    // O3: s = 46
                    if (s == 46 or s == 47 or s == 51) 
                      output_file_date << street->GetStreetConcentration(s) << "\t" ;
                  }
              }
            output_file_date << endl;
          }
        output_file_date.close();
      }
  }

  //! Initialize the output saver.
  template<class T>
  void StreetNetworkTransport<T>::InitOutputSaver()
  {
    /*** Output configuration ***/

    this->config.SetSection("[output]");
    this->config.PeekValue("Configuration_file", output_config);
    ConfigStream output_stream(output_config);
    output_stream.PeekValue("Result_dir", output_dir);
    // Section "[save]".
    output_stream.SetSection("[save]");
    output_stream.PeekValue("Text_file", text_file);
  }

  //! Sets the time step.
  /*!
    \param delta_t time step in seconds.
  */
  template<class T>
  void StreetNetworkTransport<T>::SetDelta_t(T delta_t)
  {
    this->Delta_t = delta_t;
  }

  //! Returns the time step.
  /*!
  */
  template<class T>
  T StreetNetworkTransport<T>::GetDelta_t() const
  {
    return this->Delta_t;
  }

  //! Sets the start date.
  /*!
    \param date_min start date.
  */
  template<class T>
  void StreetNetworkTransport<T>::SetDateMin(Date date_min)
  {
    this->Date_min = date_min;
  }

  //! Returns the start date.
  /*!
  */
  template<class T>
  Date StreetNetworkTransport<T>::GetDateMin() const
  {
    return this->Date_min;
  }

  // //LL: Remove stationary regime

  //! This subroutine solves a system of Ordinary Differential Equations  
  //!  with the Explicit Trapezoidal Rule algorithm (ETR).
  /*!
    \param transfer_velocity Transfer velocity between the roof and the atmosphere [m/s]
    \param temp Transfer velocity x Width x length [m3/s]
    \param outgoing_flux Outgoing flux in a street [m3/s]
    \param street_volume Street volume (lenght x width x hight) [m3]
    \param init_street_conc Street concentration at the previous time step [mug/m3]
    \param street_conc Street concentration of a determined species [mug/m3]
    \param emission_rate Emission rate of a determined species in a street [mug/s]
    \param inflow_rate Inflow rate in a street [mug/s]
    \param deposition_flux Deposition flux in a street [m3/s]
    \param conc_bg Background concentration of a determined species in a street [mug/m3]
    \param street_conc_new New street concentration of a determined species [mug/m3]
  */
  template<class T>
  void StreetNetworkTransport<T>
  ::ETRConcentration(const T transfer_velocity,
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
		     const Array<T, 1> initial_deposited_mass,
		     Array<T, 1>& new_concentration_array,
		     const T sub_delta_t,
		     const int Nsp,
		     const int bin,
		     Array<T, 1>& new_street_surface_deposited_mass)
  {
    T dtetr;
    Array<T, 1> dc1dt(Nsp);
    Array<T, 1> dc2dt(Nsp);
    dc1dt = 0.0;
    dc2dt = 0.0;
    Array<T, 1> dm1dt(Nsp);
    Array<T, 1> dm2dt(Nsp);
    Array<T, 1> street_surface_deposited_mass_tmp(Nsp);
    dm1dt = 0.0;
    dm2dt = 0.0;
    street_surface_deposited_mass_tmp = 0.0;
   
    //First time step
    CalculDCDT(transfer_velocity,
	       temp,
	       outgoing_flux,
	       street_volume,
	       concentration_array,
	       background_concentration_array,
	       emission_rate_array,
	       inflow_rate_array,
	       deposition_flux_array,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       initial_deposited_mass,
	       dc1dt,
	       Nsp,
	       bin);

    CalculDMDT(concentration_array,
	       initial_deposited_mass,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       dm1dt,
	       Nsp,
	       bin);
    
    //Update concentrations
    for (int s = 0; s < Nsp; s++)
      {
	concentration_array_tmp(s) = concentration_array(s) + dc1dt(s)*sub_delta_t;
	concentration_array_tmp(s) = max(concentration_array_tmp(s), 0.0);
	street_surface_deposited_mass_tmp(s) = initial_deposited_mass(s) + dm1dt(s)*sub_delta_t;
	street_surface_deposited_mass_tmp(s) = max(street_surface_deposited_mass_tmp(s), 0.0);
      }

    //Second time step
    CalculDCDT(transfer_velocity,
	       temp,
	       outgoing_flux,
	       street_volume,
	       concentration_array_tmp,
	       background_concentration_array,
	       emission_rate_array,
	       inflow_rate_array,
	       deposition_flux_array,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       street_surface_deposited_mass_tmp,
	       dc2dt,
	       Nsp,
	       bin);

    CalculDMDT(concentration_array_tmp,
	       street_surface_deposited_mass_tmp,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       dm2dt,
	       Nsp,
	       bin);
    
    dtetr = sub_delta_t/2.0;
    for (int s = 0; s < Nsp; s++)
      {
	new_concentration_array(s) = concentration_array(s) + (dc1dt(s) + dc2dt(s))*dtetr;
	new_concentration_array(s) = max(new_concentration_array(s), 0.0);
	new_street_surface_deposited_mass(s) = initial_deposited_mass(s) + (dm1dt(s) + dm2dt(s))*dtetr;
	new_street_surface_deposited_mass(s) = max(new_street_surface_deposited_mass(s), 0.0);
      }
  }

  //! This subroutine solves a system of Ordinary Differential Equations  
  //!  with the Explicit Trapezoidal Rule algorithm (ETR).
  /*!
    \param transfer_velocity Transfer velocity between the roof and the atmosphere [m/s]
    \param temp Transfer velocity x Width x length [m3/s]
    \param outgoing_flux Outgoing flux in a street [m3/s]
    \param street_volume Street volume (lenght x width x hight) [m3]
    \param init_street_conc Street concentration at the previous time step [mug/m3]
    \param street_conc Street concentration of a determined species [mug/m3]
    \param emission_rate Emission rate of a determined species in a street [mug/s]
    \param inflow_rate Inflow rate in a street [mug/s]
    \param deposition_flux Deposition flux in a street [m3/s]
    \param conc_bg Background concentration of a determined species in a street [mug/m3]
    \param street_conc_new New street concentration of a determined species [mug/m3]
  */
  template<class T>
  void StreetNetworkTransport<T>
  ::CalculDCDT(const T transfer_velocity,
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
	       Array<T, 1>& dcdt,
	       const int Nsp,
	       const int bin)
  {
    if (Nsp != this->Ns_aer)//gas phase and number
      for (int s = 0; s < Nsp; s++)
	{
	  dcdt(s) = (emission_rate_array(s) + inflow_rate_array(s) + resuspension_factor * street_deposited_mass_array(s) - washoff_factor_array(s) * street_deposited_mass_array(s) - outgoing_flux*concentration_array(s) - deposition_flux_array(s)*concentration_array(s) - scavenging_flux_array(s)*concentration_array(s) - temp*(concentration_array(s) - background_concentration_array(s)))/street_volume;
	}
    else //aerosol mass
      for (int s = 0; s < Nsp; s++)
	{
	  dcdt(s) = (emission_rate_array(s) + inflow_rate_array(s) + resuspension_factor * street_deposited_mass_array(s) - washoff_factor_array(s) * street_deposited_mass_array(s) - outgoing_flux*concentration_array(s) - deposition_flux_array(bin)*concentration_array(s) - scavenging_flux_array(bin)*concentration_array(s) - temp*(concentration_array(s) - background_concentration_array(s)))/street_volume;
	}
  }

  //    */
  template<class T>
  void StreetNetworkTransport<T>
  ::CalculDMDT(const Array<T, 1> concentration_array,
	       const Array<T, 1> deposited_mass_array,
	       const Array<T, 1> street_deposition_flux_array,
	       const Array<T, 1> street_scavenging_flux_array,
	       const T resuspension_factor,
	       const Array<T, 1> washoff_factor_array,
	       Array<T, 1>& dmdt,
	       const int Nsp,
	       const int bin)
  {
    if (Nsp != this->Ns_aer)// number	
      for (int s = 0; s < Nsp; s++)
	dmdt(s) = (street_deposition_flux_array(s) + street_scavenging_flux_array(s)) * concentration_array(s) - resuspension_factor * deposited_mass_array(s) - washoff_factor_array(s) * deposited_mass_array(s);
      
    else //aerosol-phase
      for (int s = 0; s < Nsp; s++)
	dmdt(s) = (street_deposition_flux_array(bin) + street_scavenging_flux_array(bin)) * concentration_array(s) - resuspension_factor * deposited_mass_array(s) - washoff_factor_array(s) * deposited_mass_array(s);
	  
  }

  //! This subroutine solves a system of Ordinary Differential Equations  
  //!  with the Explicit Trapezoidal Rule algorithm (ETR).
  /*!
    \param transfer_velocity Transfer velocity between the roof and the atmosphere [m/s]
    \param temp Transfer velocity x Width x length [m3/s]
    \param outgoing_flux Outgoing flux in a street [m3/s]
    \param street_volume Street volume (lenght x width x hight) [m3]
    \param init_street_conc Street concentration at the previous time step [mug/m3]
    \param street_conc Street concentration of a determined species [mug/m3]
    \param emission_rate Emission rate of a determined species in a street [mug/s]
    \param inflow_rate Inflow rate in a street [mug/s]
    \param deposition_flux Deposition flux in a street [m3/s]
    \param conc_bg Background concentration of a determined species in a street [mug/m3]
    \param street_conc_new New street concentration of a determined species [mug/m3]
  */
  template<class T>
  void StreetNetworkTransport<T>
  ::InitStep(T& sub_delta_t,
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
	     const int bin)
  {
    Array<T, 1> dcdt(Nsp);
    dcdt = 0.0;
    Array<T, 1> sub_delta_t_sp(Nsp);
    sub_delta_t_sp = 0.0;
    Array<T, 1> concentration_array_tmp(Nsp);
    concentration_array_tmp = 0.0;
    Array<T, 1> deposition_flux_tmp(Nsp);
    deposition_flux_tmp = 0.0;
    
    if(this->current_date == this->Date_min)
      {
	Array<T, 1> dcdt0(Nsp);
	dcdt0 = 0.0;
	Array<T, 1> dcdtmin(Nsp);
	dcdtmin = 0.0;

	CalculDCDT(transfer_velocity,
		   temp,
		   outgoing_flux,
		   street_volume,
		   concentration_array,
		   background_concentration_array,
		   emission_rate_array,
		   inflow_rate_array,
		   deposition_flux_array,
		   street_deposition_flux_array,
		   scavenging_flux_array,
		   resuspension_factor,
		   washoff_factor_array,
		   street_deposited_mass_array,
		   dcdt0,
		   Nsp,
		   bin);
	
	for (int s = 0; s < Nsp; s++)
	  {
	    concentration_array_tmp(s) = concentration_array(s) + dcdt0(s)*sub_delta_t_min;
	    concentration_array_tmp(s) = max(0.0, concentration_array_tmp(s));
	  }
	
	CalculDCDT(transfer_velocity,
		   temp,
		   outgoing_flux,
		   street_volume,
		   concentration_array_tmp,
		   background_concentration_array,
		   emission_rate_array,
		   inflow_rate_array,
		   deposition_flux_array,
		   street_deposition_flux_array,
		   scavenging_flux_array,
		   resuspension_factor,
		   washoff_factor_array,
		   street_deposited_mass_array,
		   dcdtmin,
		   Nsp,
		   bin);
	
	for (int s = 0; s < Nsp; s++)
          {
            if (dcdtmin(s) == 0.0)
              sub_delta_t_sp(s) = sub_delta_t_min;
            else
              sub_delta_t_sp(s) = abs(dcdt0(s)/dcdtmin(s) * sub_delta_t_min);
          }
      }
    else
      {
	CalculDCDT(transfer_velocity,
		   temp,
		   outgoing_flux,
		   street_volume,
		   concentration_array,
		   background_concentration_array,
		   emission_rate_array,
		   inflow_rate_array,
		   deposition_flux_array,
		   street_deposition_flux_array,
		   scavenging_flux_array,
		   resuspension_factor,
		   washoff_factor_array,
		   street_deposited_mass_array,
		   dcdt,
		   Nsp,
		   bin);
	
	for (int s = 0; s < Nsp; s++)
	  {
	    T tmp = concentration_array(s) * dcdt(s);
	    if (tmp != 0.0)
	      {
		T tscale = 0.1 * concentration_array(s)/abs(dcdt(s));
		sub_delta_t_sp(s) = min(tscale, this->Delta_t);
		sub_delta_t_sp(s) = max(sub_delta_t_sp(s), sub_delta_t_min);
	      }
	  }
      }
    sub_delta_t = min(sub_delta_t_sp);
    sub_delta_t = max(sub_delta_t, sub_delta_t_min);
    sub_delta_t = min(sub_delta_t, this->Delta_t);
  }

  template<class T>
  void StreetNetworkTransport<T>
  ::AdaptTimeStep(const Array<T, 1> new_concentration_array,
		  const Array<T, 1>concentration_array_tmp,
		  const T sub_delta_t_init,
		  const T sub_delta_t_min,
		  const T sub_delta_t_max,
		  T& sub_delta_t)
  {
    T tmp, R;
    T EPSER = 0.01; //relative error precision
    //zero init
    T n2err = 0.0;
    
    //local error estimation
    for (int s = 0; s < this->Ns; s++)
      if(new_concentration_array(s) > 0.0)
	{
	  tmp = (new_concentration_array(s) - concentration_array_tmp(s))/new_concentration_array(s);
	  n2err = n2err + tmp*tmp;
	}
    n2err = sqrt(n2err);

    //******compute new time step
    // ! first we constrain norm2 error
    // ! in order to prevent division by zero
    // ! and to keep new time step between
    // ! sub_delta_t_min and Delta_t defined 

    R = (1.0e2/1.0e-5);
    tmp = R*R;
    n2err = min(n2err, EPSER*tmp);
    n2err = max(EPSER/tmp, n2err);

    //formula to compute new time step
    sub_delta_t = sub_delta_t_init*sqrt(EPSER/n2err);
    sub_delta_t = min(sub_delta_t, sub_delta_t_max);
    sub_delta_t = max(sub_delta_t, sub_delta_t_min);
  }

  //! Sets the dry deposition flux for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetStreetDryDepositionFlux()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
	  for (int j = 0; j < Ns_dep; ++j)
	    if (this->species_list[s] == species_list_dep[j])
	      {
		StreetDryDepositionFlux(s, ist) = street->GetStreetDryDepositionFlux(s);
	      }
        ++ist;
      }
  }

/*
//! Sets the wall dry deposition flux for all species at the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetWallDryDepositionFlux()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
	  for (int j = 0; j < Ns_dep; ++j)
	    if (this->species_list[s] == species_list_dep[j])
	      {
		WallDryDepositionFlux(s, ist) = street->GetWallDryDepositionFlux(s);
	      }
        ++ist;
      }
  }
*/
  //! Sets the dry deposition rate for the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetStreetDryDepositionRate()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
	  for (int j = 0; j < Ns_dep; ++j)
	    if (this->species_list[s] == species_list_dep[j])
	      {
		StreetDryDepositionRate(s, ist) = street->GetStreetDryDepositionFlux(s) * street->GetLength() * street->GetWidth();
	      }
        ++ist;
      }
  }

  //! Sets the wall dry deposition rate for all species at the whole street-network
  template<class T>
  void StreetNetworkTransport<T>::SetWallDryDepositionRate()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
	  for (int j = 0; j < Ns_dep; ++j)
	    if (this->species_list[s] == species_list_dep[j])
	      {
		WallDryDepositionRate(s, ist) = street->GetWallDryDepositionFlux(s) *
		  street->GetHeight() * street->GetLength() * 2.0;
	      }
        ++ist;
      }
  }  

  //! Sets the scavening flux for all species in the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetStreetScavengingFlux()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
	  for (int j = 0; j < Ns_scav; ++j)
	    if (this->species_list[s] == species_list_scav[j])
	      {
		StreetScavengingFlux(s, ist) = street->GetStreetScavengingFlux(s);
	      }
        ++ist;
      }
  }

  //! Sets the scavening rate for all species in the whole street-network.
  template<class T>
  void StreetNetworkTransport<T>::SetStreetScavengingRate()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns; ++s)
	  for (int j = 0; j < Ns_scav; ++j)
	    if (this->species_list[s] == species_list_scav[j])
	      {
		StreetScavengingRate(s, ist) = street->GetStreetScavengingFlux(s) *
		  street->GetWidth() * street->GetLength();
	      }
        ++ist;
      }
  }

  // Rosenbrock method
  template<class T>
  void StreetNetworkTransport<T>
  ::RosenbrockConcentration(const T transfer_velocity,
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
			    const int bin)
  {
    Array<T, 1> dc1dt(Nsp); // Differential term dc/dt at initial time (sub time step)
    Array<T, 1> dc2dt(Nsp); // Differential term dc/dt at final time (sub time step)
    Array<T, 1> J(Nsp);
    Array<T, 1> k1(Nsp);
    Array<T, 1> k2(Nsp);

    
    Array<T, 1> dm1dt(Nsp); // Differential term dm/dt at initial time (sub time step)
    Array<T, 1> dm2dt(Nsp); // Differential term dm/dt at final time (sub time step)
    Array<T, 1> J_surf(Nsp);
    Array<T, 1> k1_surf(Nsp);
    Array<T, 1> k2_surf(Nsp);
    Array<T, 1> street_surface_deposited_mass_tmp(Nsp);
    
    
    T gamma;
    
    dc1dt = 0.0;
    dc2dt = 0.0;
    J = 0.0;
    k1 = 0.0;
    k2 = 0.0;
    
    dm1dt = 0.0;
    dm2dt = 0.0;
    J_surf = 0.0;
    k1_surf = 0.0;
    k2_surf = 0.0;
    street_surface_deposited_mass_tmp = 0.0;
    
    gamma = 1 + 1 / sqrt(2);
    
    // Compute the differential term dc/dt at initial time
    CalculDCDT(transfer_velocity,
	       temp,
	       outgoing_flux,
	       street_volume,
	       concentration_array,
	       background_concentration_array,
	       emission_rate_array,
	       inflow_rate_array,
	       deposition_flux_array,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       initial_deposited_mass,
	       dc1dt,
	       Nsp,
	       bin);
    
    // Compute the Jacobian at initial time - mass balance in the street volume
    Jacobian(inflow_flux,
	     outgoing_flux,
	     temp,
	     deposition_flux_array,
	     J,
	     street_volume,
	     Nsp,
	     bin);

    // Compute mass variation on the surface
    CalculDMDT(concentration_array,
	       initial_deposited_mass,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       dm1dt,
	       Nsp,
	       bin);
    // Compute the Jacobian at initial time - mass balance in the street surface
    Jacobian(resuspension_factor,
	     washoff_factor_array,
	     deposition_flux_array,
	     J_surf,
	     street_volume,
	     Nsp,
	     bin);    

    // Compute k1 - street volume
    for (int s = 0; s < Nsp; s++)
      k1(s) = dc1dt(s) / (1 - gamma * sub_delta_t * J(s));

    // Compute temporary concentrations at final time
    for (int s = 0; s < Nsp; s++)
      {
	concentration_array_tmp(s) = concentration_array(s) + k1(s) * sub_delta_t;
	concentration_array_tmp(s) = max(concentration_array_tmp(s), 0.0);
      }
    // Compute k1 - street surface
    for (int s = 0; s < Nsp; s++)
      k1_surf(s) = dm1dt(s) / (1 - gamma * sub_delta_t * J_surf(s));

    // Compute temporary concentrations at final time
    for (int s = 0; s < Nsp; s++)
      {
	street_surface_deposited_mass_tmp(s) = initial_deposited_mass(s) + k1_surf(s) * sub_delta_t;
	street_surface_deposited_mass_tmp(s) = max(street_surface_deposited_mass_tmp(s), 0.0);
      }
  
    // Compute the differential term dc/dt at final time
    CalculDCDT(transfer_velocity,
	       temp,
	       outgoing_flux,
	       street_volume,
	       concentration_array_tmp,
	       background_concentration_array,
	       emission_rate_array,
	       inflow_rate_array,
	       deposition_flux_array,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       street_surface_deposited_mass_tmp,
	       dc2dt,
	       Nsp,
	       bin);

    CalculDMDT(concentration_array_tmp,
	       street_surface_deposited_mass_tmp,
	       street_deposition_flux_array,
	       scavenging_flux_array,
	       resuspension_factor,
	       washoff_factor_array,
	       dm2dt,
	       Nsp,
	       bin);
    
    // Compute k2
    for (int s = 0; s < Nsp; s++)
      k2(s) = (dc2dt(s) - 2 * k1(s)) / (1 - gamma * sub_delta_t * J(s));
      
    // Compute concentrations at final time
    for (int s = 0; s < Nsp; s++)
      {
	new_concentration_array(s) = concentration_array(s) + (3 * k1(s) + k2(s)) * sub_delta_t / 2;
	new_concentration_array(s) = max(new_concentration_array(s), 0.0);
      }
    
    // Compute k2
    for (int s = 0; s < Nsp; s++)
      k2_surf(s) = (dm2dt(s) - 2 * k1_surf(s)) / (1 - gamma * sub_delta_t * J_surf(s));
      
  }
  
  // Compute Jacobian matrix (1x1) by Euler
  template<class T>
  void StreetNetworkTransport<T>
  ::Jacobian(const T inflow_flux,
	     const T outgoing_flux,
	     const T temp,
	     const Array<T, 1> deposition_flux_array,
	     Array<T, 1>& J,
	     const T street_volume,
	     int Nsp,
	     int bin)
  {
    for (int s = 0; s < Nsp; s++)
      {
	if(Nsp != this->Ns_aer) //gas-phase pollutants and number
	  J(s) = (inflow_flux - temp - outgoing_flux - deposition_flux_array(s)) / street_volume;
	else //particle-phase pollutants
	  J(s) = (inflow_flux - temp - outgoing_flux - deposition_flux_array(bin)) / street_volume; 
      }
  }

  // Compute Jacobian matrix (1x1) by Euler
  template<class T>
  void StreetNetworkTransport<T>
  ::Jacobian(const T resuspension_factor,
	     const Array<T, 1> washoff_factor_array,
	     const Array<T, 1> deposition_flux_array,
	     Array<T, 1>& J,
	     const T street_volume,
	     int Nsp,
	     int bin)
  {
    for (int s = 0; s < Nsp; s++)
      {
	if (Nsp != this->Ns_aer) //gas-phase pollutants and number
	  J(s) = deposition_flux_array(s)/street_volume - resuspension_factor - washoff_factor_array(s);
	else //particle-phase
	  J(s) = deposition_flux_array(bin)/street_volume - resuspension_factor - washoff_factor_array(s);
      }
  }
  

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKTRANSPORT_CXX
#endif
