#ifndef POLYPHEMUS_FILE_MODELS_STREETNETWORKAEROSOL_CXX

//////////////
// INCLUDES //

#include "StreetNetworkAerosol.hxx"
#include <vector>
#include <iostream>
#include <numeric>


// INCLUDES //
//////////////
  
namespace Polyphemus
{

  template<class T, class ClassChemistry>
  const T StreetNetworkAerosol<T, ClassChemistry>::pi = acos(-1);

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds the model. Nothing else is performed.
  */
  template<class T, class ClassChemistry>
  StreetNetworkAerosol<T, ClassChemistry>::StreetNetworkAerosol():
  StreetNetworkChemistry<T, ClassChemistry>()
  {
  }


  //! Main constructor.
  /*!
    \param config_file configuration filename.
  */
  template<class T, class ClassChemistry>
  StreetNetworkAerosol<T, ClassChemistry>::StreetNetworkAerosol(string config_file):
    StreetNetworkChemistry<T, ClassChemistry>(config_file)
  {
    this->D3_map["StreetSurfaceDepositedMass_aer"] = &StreetSurfaceDepositedMass_aer;
    this->D3_map["StreetResuspensionRate_aer"] = &StreetResuspensionRate_aer;
    this->D3_map["StreetSurfaceDryDepositionRate_aer"] = &StreetSurfaceDryDepositionRate_aer;
    this->D3_map["StreetWashoffRate_aer"] = &StreetWashoffRate_aer;
  }
  

  //! Destructor.
  template<class T, class ClassChemistry>
  StreetNetworkAerosol<T, ClassChemistry>::~StreetNetworkAerosol()
  {
    this->ClearStreetVector();
    this->ClearIntersectionVector();
  }


  //! Append a bin number in the bin vector.
  /*!
    \param bin_list_aer bin vector.
    \param b bin number.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>
  ::AppendBinList_aer(int b, vector<int>& bin_list_aer)
  {

    if (bin_list_aer.empty())
      bin_list_aer.push_back(b);
    else
      {
	vector<int>::iterator pos;
		
	pos = find(bin_list_aer.begin(),
                   bin_list_aer.end(), b);
        
	if (pos == bin_list_aer.end())	
	  {
            int ins = 0;
            for (int b_id = 0; b_id < int(bin_list_aer.size()); ++b_id)
              if (b > bin_list_aer[b_id])
                ++ins;
	    bin_list_aer.insert(bin_list_aer.begin() + ins, b);
	  }
      }  
  }
  
  
  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::ReadConfiguration()
  {
    StreetNetworkChemistry<T, ClassChemistry>::ReadConfiguration();
    /*** Options ***/
    
    this->config.SetSection("[options]");

    this->config.PeekValue("With_initial_condition_aerosol",
			   this->option_process
			   ["with_initial_condition_aer"]);
    this->config.PeekValue("With_initial_condition_number_aerosol",
			   this->option_process
			   ["with_initial_condition_number_aer"]);
    
    //! Deposition
    this->config.PeekValue("With_deposition_aerosol",
			   this->option_process["with_deposition_aer"]);
    this->config.PeekValue("With_drainage_aerosol",
			   this->option_process["with_drainage_aer"]);
    this->config.PeekValue("Min_water_drainage",
			   "> 0.0",
			   min_water_drainage); //mm
    
    if (this->option_process["with_deposition_aer"])
      {
        if (this->config.Check("Collect_dry_flux_aerosol"))
          this->config.PeekValue("Collect_dry_flux_aerosol",
                                 this->option_process["collect_dry_flux_aer"]);
        else
          this->option_process["collect_dry_flux_aer"] = false;
        this->config.PeekValue("Particles_dry_velocity_option",
			       "zhang|giardina",
			       particles_dry_velocity_option);
        this->config.PeekValue("Brownian_diffusion_resistence_option",
                                 "zhang|giardina|seigneur|chamberlain",
                                 brownian_diffusion_resistence_option);
      }
    
    //! Scavenging
    this->config.PeekValue("With_scavenging_aerosol",
			   this->option_process["with_scavenging_aer"]);
    
    if (this->option_process["with_scavenging_aer"])
      {
        if (this->config.Check("Collect_wet_flux_aerosol"))
          this->config.PeekValue("Collect_wet_flux_aerosol",
                                 this->option_process["collect_wet_flux_aer"]);
        else
          this->option_process["collect_wet_flux_aer"] = false;
      }
    
    this->config.PeekValue("With_resuspension",
                           this->option_process["with_resuspension"]);

    //! Resuspension
    if (this->option_process["with_resuspension"])
      {

        if (this->config.Check("Resuspension"))
          this->config.PeekValue("Resuspension", file_resuspension);
        else
          file_resuspension = "resuspension.dat";
        ConfigStream config_resuspension(file_resuspension);
        //Precipitation that indicates a complete drainage in streets
        this->config.PeekValue("Max_rain", max_rain);
	
        config_resuspension.PeekValue("f0_hdv", 
                                      f0_hdv);
        config_resuspension.PeekValue("f0_lcv",
                                      f0_lcv);
        config_resuspension.PeekValue("f0_pc",
                                      f0_pc);
        config_resuspension.PeekValue("f0_2R",
                                      f0_2R);
        config_resuspension.PeekValue("mean_speed_2R",
                                      mean_speed_2R);
        config_resuspension.PeekValue("mean_speed_HDV",
                                      mean_speed_HDV);
        config_resuspension.PeekValue("mean_speed_PC",
                                      mean_speed_PC);
        config_resuspension.PeekValue("mean_speed_LCV",
                                      mean_speed_LCV);
        config_resuspension.PeekValue("mean_speed_highway_2R",
                                      mean_speed_highway_2R);
        config_resuspension.PeekValue("mean_speed_highway_HDV",
                                      mean_speed_highway_HDV);
        config_resuspension.PeekValue("mean_speed_highway_PC",
                                      mean_speed_highway_PC);
        config_resuspension.PeekValue("mean_speed_highway_LCV",
                                      mean_speed_highway_LCV);
      }
   
    
    //! Number concentration
    this->config.PeekValue("With_number_concentration",
                           this->option_process["with_number_concentration"]);
    if (this->option_process["with_number_concentration"])
      {
        if (this->config.Check("Number_computation_option"))
          this->config.PeekValue("Number_computation_option",
                                 "based_on_mass|based_on_transport",
                                 number_computation_option);
        else
          number_computation_option = "based_on_mass"; // recommended

        if (this->config.Check("With_bg_number_concentration_data"))        
          this->config.PeekValue("With_bg_number_concentration_data",
                                 this->option_process["with_bg_number_concentration_data"]);
        else
          this->option_process["with_bg_number_concentration_data"] = false;
    
        if (this->config.Check("With_emission_number_data"))
          this->config.PeekValue("With_emission_number_data",
                                 this->option_process["with_emission_number_data"]);
        else
          this->option_process["with_emission_number_data"] = false;

      }
    
    //! Density in kg / m^3.
    this->config.PeekValue("With_fixed_density",
			   this->option_process["with_fixed_density"]);    
    this->config.PeekValue("Fixed_aerosol_density", "positive",
			   fixed_density_aer);

    //! Option if wet diameter is calculated using Gerber or in chemistry
    this->config.PeekValue("Wet_diameter_option",
			   "gerber|chemistry|none",
			   wet_diameter_option); 

    
    /* Options from ssh-aerosol */
    
    this->Chemistry_.InitSharedLib(*this);
        
    ssh_Ns_aer_nolayer = this->Chemistry_.GetNsAerNoLayer();        
    this->Ns_aer = this->Chemistry_.GetNsAer();
    N_layer = this->Chemistry_.GetNLayer();
    i_hydrophilic = this->Chemistry_.GetIHydrophilic();        
    ssh_mass_density_layers.resize(this->Ns_aer);
    // conversion to g/cm3
    ssh_mass_density_layers = this->Chemistry_.GetMassDensityLayers() * 1.0e6;

    index_species_layers.resize(this->Ns_aer);

    this->Nbin_aer = this->Chemistry_.GetNbinAer();
    this->option_process["with_external_composition"] =
      this->Chemistry_.IsExternalComposition();
        
    this->Nsize_section_aer = this->Chemistry_.GetNSizeSectionAer();

    //! Reads mass density aerosol
    Rho_species_aer.resize(this->Ns_aer);
    for (int i = 0; i < this->Ns_aer; i++)
      Rho_species_aer[i] = ssh_mass_density_layers(i);
   
    /*** Species and bins ***/

    this->config.SetSection("[domain]");
    // Reads bin bounds.
    this->config.Find("Bin_bounds");
    bin_list = split(this->config.GetLine());
    this->Nsize_section_aer = int(bin_list.size()) - 1;
    this->config.PeekValue("Species", this->file_species);
    // Opens the file that describes species.
    ConfigStream species_stream(this->file_species);
    vector<string> species_bin;
    string species;
    int bin_index;

    
    if(this->option_process["with_drainage_aer"])
      {
	species_stream.SetSection("[drainage_efficiency]");
	while (!species_stream.IsEmpty())
	  {
	    species = species_stream.GetElement();
	    species_stream.GetNumber(drainage_efficiency[species]);
	  }
      }
    
    //! Get aerosol types which are defined in SSH-aerosol.
    aerosol_type.resize(ssh_Ns_aer_nolayer);
    aerosol_type = this->Chemistry_.GetAerosolType();

    //! Modify the number of aerosol species
    //! considering additional phases from SSH-aerosol.
    /*
      \param species_list_aer_nolayer species name
      \param species_list_aer species name + postfix
     */
    
    n_nonorganic = 0;
    n_organics = 0;
    int Ns_aer_check = 0;
    int total_layer = N_layer + i_hydrophilic;

    //! Get aerosol types which are defined in SSH-aerosol.
    species_list_aer_nolayer.resize(ssh_Ns_aer_nolayer);
    species_list_aer_nolayer = this->Chemistry_.GetAerosolSpecName();
    
    for (int s = 0; s < ssh_Ns_aer_nolayer; s++)
      {
        int index = 0;
        string species_name = species_list_aer_nolayer(s);
        //! aerosol_type is 4 for organic aerosols in SSH-aerosol.
        if ( aerosol_type(s) != 4)
          {
            this->species_list_aer.push_back(species_name);
            ++Ns_aer_check;
            index_species_layers(index) = index;
            ++index;
            ++n_nonorganic;
          }
        else
          {
            /* A postfix "_layer(+ number)" is added in the species names 
               in the layers. The species without the postfix correspond to 
               the species in the first layer. Input data given in 
               the initial condition and boundary condition are used for 
               these species without the postfix.
             */
            for (int l = 0; l < total_layer; l++)
              {
                string species_name_layer;
                if (i_hydrophilic == 1 and l == total_layer - 1)
                  species_name_layer = species_name + "_hydrophilic";
                else
                  species_name_layer = species_name + "_layer" + to_str(l + 1);
                ++n_organics;
                ++Ns_aer_check;
                index_species_layers(index + l) = index;
                if (l == 0)
                  this->species_list_aer.push_back(species_name);
                else
                  this->species_list_aer.push_back(species_name_layer);
              }
            ++index;
          }
      }

    if (Ns_aer_check != this->Ns_aer or
        Ns_aer_check != (n_nonorganic + n_organics))
      throw string("Something is wrong in the list of aerosol species.") +
        to_str(Ns_aer_check) + " "  + to_str(this->Ns_aer);
    
    
    this->list_aer = this->species_list_aer;
    // Add "Number" to the species list aerosol
    if (this->option_process["with_number_concentration"])
      this->list_aer.push_back("Number");


    aerosol_species_group_relation.resize(this->Ns_aer);

    this->Ncomposition_aer = this->Chemistry_.GetNcomposition_aer();
      
    this->Ngroup_aer = this->Chemistry_.GetNgroup_aer();

    composition_bounds.resize(this->Ncomposition_aer,this->Ngroup_aer,2);

    composition_bounds = this->Chemistry_.GetCompositionBounds();

    Nfraction_aer = this->Chemistry_.GetNfraction_aer();

    
    if (this->option_process["with_external_composition"])
      {

        for (int i = 0; i < this->Ns_aer; i++)
          // Conversion Fortran index n to C++ index n-1. 
          this->aerosol_species_group_relation(i) = this->Chemistry_.GetIndexGroupsExt()(i) - 1;

        //! For water
        this->aerosol_species_group_relation(this->Ns_aer - 1) = 0;        
      }
    else//in case of internal mixing
      {
        this->groups_list_aer.push_back("None");
        this->Ngroup_aer=int(this->groups_list_aer.size());
        for (int i = 0; i < this->Ns_aer; i++)
          this->aerosol_species_group_relation(i)=0;
      }

    
     /*** Input files ***/

    /*! The configuration-file path is the field "Data_description" in the main
       configuration file.     
    */
    
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    ConfigStream data_description_stream(data_description_file);

    //! Dissolution heat (kcal / mol at 298K).
    species_stream.SetSection("[dissolution_heat]");
    while (!species_stream.IsEmpty())
      {
	species = species_stream.GetElement();
	species_stream.GetNumber(dissolution_heat[species]);
      }
    
   
    //! The configuration-file path is the field "Data_description" in the main
    //! configuration file.
    this->config.SetSection("[data]");


    //! Aerosol emission files.
    this->input_files["traffic"].Read(data_description_file, "traffic");
    data_description_stream.SetSection("[traffic]");
    data_description_stream.PeekValue("Nt", "> 0", Nt_traffic);
    
    this->input_files["emission_aer"].Read(data_description_file,
                                           "emission_aer");
    data_description_stream.SetSection("[emission_aer]");
    if (data_description_stream.Check("Format"))
      data_description_stream.PeekValue("Format", aerosol_emission_format);
    else
      aerosol_emission_format = "Internal";
    data_description_stream.PeekValue("Nt", "> 0", Nt_emis_aer);
    this->input_files["emission_aer"]
      .ReadFiles(data_description_file, "emission_aer");
    
    for (map<string, string>::iterator i
           = this->input_files["emission_aer"].Begin();
         i != this->input_files["emission_aer"].End(); i++)
      {
	species_bin = split(i->first, "_");
	
	if (species_bin.size() != 2)
	  throw string("Species \"") + i->first + "\" is badly formatted.";
	species = species_bin[0];
	bin_index = convert<int>(species_bin[1]);

        //! Built bin list with initial conditions
	if (bin_index < this->Nsize_section_aer)
          //! Build with size section
          AppendBinList_aer(bin_index, emis_bin_list_aer);
	else
	  throw string("aerosol_emission: Bin index ") +
            "in data file are out of range";
	unsigned int j = 0;
	while (j < species_list_emis_aer.size()
	       && species_list_emis_aer[j].first != species)
	  j++;
	if (j == species_list_emis_aer.size())
	  species_list_emis_aer
	    .push_back(pair<string, vector<int> >(species, vector<int>()));
	species_list_emis_aer[j].second.push_back(bin_index);
      }
    Ns_emis_aer = int(species_list_emis_aer.size()); 
    Nb_emis_aer = int(emis_bin_list_aer.size());

    
    //! Background concentrations files.   
    if (this->option_process["with_local_data"])
      {
	// Aerosol species
	this->input_files["bg_concentration_aer"].Read(data_description_file,
                                                       "bg_concentration_aer");
	data_description_stream.SetSection("[bg_concentration_aer]");
        if (data_description_stream.Check("Format"))
          data_description_stream.PeekValue("Format", aerosol_bg_format);
        else
          aerosol_bg_format = "Internal";
	data_description_stream.PeekValue("Nt", "> 0", Nt_bg_aer);
	this->input_files["bg_concentration_aer"]
	  .Read(data_description_file, "bg_concentration_aer");
	
	for (map<string, string>::iterator i 
	      = this->input_files["bg_concentration_aer"].Begin();
	    i != this->input_files["bg_concentration_aer"].End(); i++)
	  {
	    species_bin = split(i->first, "_");
	    if (species_bin.size() != 2)
	      throw string("Species \"") + i->first + "\" is badly formatted.";
	    species = species_bin[0];
	    bin_index = convert<int>(species_bin[1]);
	    if (bin_index < this->Nsize_section_aer)
	      //! Built bin list with boundary conditions
	      AppendBinList_aer(bin_index, bg_bin_list_aer);
	    else
	      throw string("background_concentration_aer: ") +
                "Bin index in data file are out of range";
	    unsigned int j = 0;
	    while (j < species_list_bg_aer.size()
		  && species_list_bg_aer[j].first != species)
	      j++;
	    if (j == species_list_bg_aer.size())
	      species_list_bg_aer
		.push_back(pair<string, vector<int> >(species, vector<int>()));
	    species_list_bg_aer[j].second.push_back(bin_index);
	  }
	Ns_bg_aer = int(species_list_bg_aer.size());
	Nb_bg_aer = int(bg_bin_list_aer.size());

	//! Aerosol background number concentration
	if (this->option_process["with_number_concentration"])
	  if (this->option_process["with_bg_number_concentration_data"])
	    this->input_files["bg_number_concentration"]
	      .ReadFiles(data_description_file, "bg_number_concentration");
      } 

    //! Aerosol scavenging;
    if (this->option_process["with_scavenging_aer"])
      this->input_files["scavenging_aer"]
	.ReadFields(data_description_file, "scavenging_aerosol");
    else
      this->input_files["scavenging_aer"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["scavenging_aer"].Begin();
	 i != this->input_files["scavenging_aer"].End(); i++)
      {
	if (!is_integer(i->first))
	  throw string("In scavenging fields, \"") + i->first
	    + "\" is not a bin index.";
	bin_index = convert<int>(i->first);
	scav_bin_list_aer.push_back(bin_index);
      }
    Nb_scav_aer = int(scav_bin_list_aer.size());
    Nbin_scav_aer = Nb_scav_aer * this->Ncomposition_aer;

    //! Deposition velocities.
    if (this->option_process["with_deposition_aer"])
      this->input_files["deposition_velocity_aer"]
	.ReadFields(data_description_file, "deposition_velocity_aerosol");
    else
      this->input_files["deposition_velocity_aer"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["deposition_velocity_aer"].Begin();
	 i != this->input_files["deposition_velocity_aer"].End(); i++)
      {
	if (!is_integer(i->first))
	  throw string("In deposition fields, \"") + i->first
	    + "\" is not a bin index.";
	bin_index = convert<int>(i->first);
	dep_bin_list_aer.push_back(bin_index);
      }
    Nb_dep_aer = int(dep_bin_list_aer.size());
    Nbin_dep_aer = Nb_dep_aer * this->Ncomposition_aer;


    //! Initial conditions.
    if (this->option_process["with_initial_condition_aer"])
    {
      data_description_stream.SetSection("[initial_condition_aerosol]");
      if (this->option_process["with_external_composition"])
        {
          if (data_description_stream.Check("Format"))
            data_description_stream.PeekValue("Format", ic_format);
          else
            ic_format = "Internal";
        }
      else
        ic_format = "Internal";
      this->input_files["initial_condition_aer"]
	.ReadFiles(data_description_file, "initial_condition_aerosol");
    }
    else
      this->input_files["initial_condition_aer"].Empty();

    //! Read species list from the configuration file.
    for (map<string, string>::iterator i
           = this->input_files["initial_condition_aer"].Begin();
         i != this->input_files["initial_condition_aer"].End(); i++)
      {
	species_bin = split(i->first, "_");
	if (species_bin.size() != 2)
	  throw string("Species \"") + i->first + "\" is badly formatted.";
	species = species_bin[0];
	bin_index = convert<int>(species_bin[1]);

        //! Built bin list with initial conditions => ic_bin_list_aer
	if (bin_index < this->Nsize_section_aer)
          {
            vector<int>::iterator pos = find(ic_bin_list_aer.begin(),
                                             ic_bin_list_aer.end(), bin_index);
            if (pos == ic_bin_list_aer.end())
              ic_bin_list_aer.push_back(bin_index);            
          }
	else
	  throw string("initial_condition_aer: ") +
                       "Bin index in data file are out of range";
        
	unsigned int j = 0;
	while (j < species_list_ic_aer.size()
	       && species_list_ic_aer[j].first != species)
	  j++;
	if (j == species_list_ic_aer.size())
	  species_list_ic_aer
	    .push_back(pair<string, vector<int> >(species, vector<int>()));
	species_list_ic_aer[j].second.push_back(bin_index);
      }
    Ns_ic_aer = int(species_list_ic_aer.size());
    Nb_ic_aer = int(ic_bin_list_aer.size());

    //! Added for number
    if (this->option_process["with_number_concentration"])
      if (this->option_process["with_initial_condition_number_aer"])
	this->input_files["initial_condition_aer"]
	  .ReadFiles(data_description_file, "initial_number_aerosol");

    CheckConfiguration();
    
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI    
    MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
#else
    this->rank = 0;
#endif
    if (this->rank == 0)
      DisplayConfiguration();
  }

  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::DisplayConfiguration()
  {
    // cout << "Number of non-organic species " << n_nonorganic << endl;    
    // cout << "Number of organic species " << n_organics << endl;
    // cout << "Number of all species " << n_nonorganic + n_organics << endl;    
    // cout << "Number of aerosols species without considering layers from ssh: " << ssh_Ns_aer_nolayer << endl;
    // cout << "Number of aerosols species considering layers from ssh: " << this->Ns_aer << endl;        
    // cout << "Number of layer from ssh: " << N_layer << endl;
    // cout << "Additional layer in the hydrophilic phase: " << i_hydrophilic << endl;
    // cout << "Number of aerosol bins: " << this->Nbin_aer << endl;
    // if (this->option_process["with_external_composition"])
    //   cout << "Externally mixed." << endl;
    // else
    //   cout << "Internally mixed." << endl;
  }
  
  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::Init()
  {
    this->SetCurrentDate(this->Date_min);
    InitStreet();
    this->InitIntersection();
    Allocate();
    this->ComputeStreetAngle();
    if (this->option_process["with_chemistry"])
      InitChemistry();

    if (this->option_process["with_tree_deposition"] and
      (!this->option_process["with_deposition_aer"]))
      throw string("Please activate the with_deposition_aer option to compute particles deposition on tree leaves.");

    if (this->option_process["with_tree_aerodynamic"] or this->option_process["with_tree_deposition"])
      this->InitTree();
    if (this->option_process["with_tree_deposition"])
      this->ReadTreeParamDeposition();

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Initialize the number of sources to parallelize.
    BaseModuleParallel::InitSource(this->GetNStreet());
    //
    BaseModuleParallel::BuildPartition_source();
#endif

  }
  
  /*! Initializes background concentrations and photolysis parameters.
    \param meteo ConfigStream instance through which background concentrations
    and photolysis rates may be read.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>
  ::InitChemistry()
  {
// #ifndef POLYPHEMUS_WITH_SSH_AEROSOL                
    // this->InitPhotolysis(this->current_date, this->PhotolysisRate_f);
// #endif
    Chemistry_.Init(*this);

  }  

  //! Streets initialization.
  /*!
   */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::InitStreet()
  {  
    Array<T, 1> init_conc(this->Ns);
    init_conc = 0.0;
    Array<T, 2> init_conc_aer(this->Ns_aer, this->Nbin_aer);
    init_conc_aer = 0.0;  
    Array<T, 1> init_number_conc(this->Nbin_aer);
    init_number_conc = 0.0;    

    for (int i = 0; i < this->total_nstreet; ++i)
      {
        Street<T>* street = 
          new StreetAerosol<T>(this->id_street(i),
			       this->begin_inter(i),
			       this->end_inter(i),
			       this->length(i),
			       this->width(i),
			       this->height(i),
			       this->typo(i),
			       this->Ns,
			       this->Nr_photolysis,
			       this->Ns_aer,
			       this->Nbin_aer);

        for (int s = 0; s < this->Ns; s++)
          street->SetStreetConcentration(init_conc(s), s);
        
        for (int s = 0; s < this->Ns_aer; s++)
	  for (int b = 0; b < this->Nbin_aer; b++)
	    street->SetStreetConcentration_aer(init_conc_aer(s,b), s, b);
	  
        for (int b = 0; b < this->Nbin_aer; b++)
            street->SetStreetNumberConcentration(init_number_conc(b), b); 	  
	  
        this->StreetVector.push_back(street);
      }
    this->current_street = this->StreetVector.begin();   
  }

  //! Checks the configuration.
  /*! 
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::CheckConfiguration()
  {
    StreetNetworkChemistry<T, ClassChemistry>::CheckConfiguration();  
    if(this->option_process["with_resuspension"] and
       !this->option_process["with_deposition_aer"])
      throw string("Please activate the aerosol deposition ") +
        "to compute particles resuspension.";
    if(this->option_process["with_deposition_aer"] and
       !this->option_process["with_deposition"])
      throw string("Please activate the with_deposition ") +
        "to compute particles deposition.";

    if (!this->option_process["with_deposition_aer"]
	&& this->option_process["collect_dry_flux_aer"])
      throw string("Aerosol dry deposition fluxes cannot be collected") +
	" without deposition.";

    if (this->option_process["collect_wet_flux_aer"]
	&& (!this->option_process["with_scavenging_aer"]))
      throw string("Aerosol wet deposition fluxes cannot be collected") +
	" without scavenging.";

  }

  
  //! Checks whether aerosols in a given bin have scavenging.
  /*!
    \param b bin number.
    \return True if the aerosols in bin \a b has deposition velocities, false
    otherwise.
  */
  template<class T, class ClassChemistry>
  bool StreetNetworkAerosol<T, ClassChemistry>
  ::HasScavenging_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    return find(scav_bin_list_aer.begin(), scav_bin_list_aer.end(), b)
      != scav_bin_list_aer.end();
  }


  //! Returns the index in scavenging of a given aerosol bin.
  /*!
    \param b bin number.
    \return The aerosol bin index in scavenging.
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>
  ::ScavengingIndex_aer(int b) const
  {
    return find(scav_bin_list_aer.begin(), scav_bin_list_aer.end(), b)
      - scav_bin_list_aer.begin();
  }

  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::Allocate()
  {
    StreetNetworkChemistry<T, ClassChemistry>::Allocate();

    GridS3D_aer = RegularGrid<T>(this->Ns_aer);
    GridB3D_aer = RegularGrid<T>(this->Nbin_aer);
    
    GridS2D_aer = RegularGrid<T>(this->Ns_aer);
    GridB2D_aer = RegularGrid<T>(this->Nbin_aer);

    GridST3D = RegularGrid<T>(this->total_nstreet);
    
    BinBound_aer.resize(this->Nsize_section_aer + 1);
    // Reads bin bounds in micrometers and converts it to meters.
    for (int i = 0; i < this->Nsize_section_aer + 1; i++)
      BinBound_aer(i) = 1.e-6 * convert<T>(bin_list[i]); 
    
    /*** Emission_aer ***/

    GridS_emis_aer = RegularGrid<T>(Ns_emis_aer);
    GridB_emis_aer = RegularGrid<T>(Nb_emis_aer*this->Ncomposition_aer);   

    /*** Initial condition ***/
    GridB2D_aer_i = RegularGrid<T>(this->Nsize_section_aer);
    GridB3D_aer_i = RegularGrid<T>(this->Nsize_section_aer);    
    
    if (this->option_process["with_local_data"])
      {
	LiquidWaterContent_i.Resize(this->GridST2D);
	LiquidWaterContent_f.Resize(this->GridST2D);
	FileLiquidWaterContent_i.Resize(this->GridST2D);
	FileLiquidWaterContent_f.Resize(this->GridST2D);
	
	LiquidWaterContent_i.SetZero();
	LiquidWaterContent_f.SetZero();
	FileLiquidWaterContent_i.SetZero();
	FileLiquidWaterContent_f.SetZero();
      }

    RelativeHumidity_i.Resize(this->GridST2D);
    RelativeHumidity_f.Resize(this->GridST2D);
    RelativeHumidity_i.SetZero();
    RelativeHumidity_f.SetZero();
    
    Emission_aer_i.Resize(GridS_emis_aer, GridB_emis_aer, this->GridST2D);
    Emission_aer_f.Resize(GridS_emis_aer, GridB_emis_aer, this->GridST2D);
    FileEmission_aer_i.Resize(GridS_emis_aer, GridB_emis_aer, this->GridST2D);
    FileEmission_aer_f.Resize(GridS_emis_aer, GridB_emis_aer, this->GridST2D);
    NumberEmission_aer_i.Resize(GridB_emis_aer, this->GridST2D);
    NumberEmission_aer_f.Resize(GridB_emis_aer, this->GridST2D);
    FileNumberEmission_aer_i.Resize(GridB_emis_aer, this->GridST2D);
    FileNumberEmission_aer_f.Resize(GridB_emis_aer, this->GridST2D);

    Emission_aer_i.SetZero();
    Emission_aer_f.SetZero();
    FileEmission_aer_i.SetZero();
    FileEmission_aer_f.SetZero();
    NumberEmission_aer_i.SetZero();
    NumberEmission_aer_f.SetZero();
    FileNumberEmission_aer_i.SetZero();
    FileNumberEmission_aer_f.SetZero();
    
    /*** Background_aer ***/
    if (this->option_process["with_local_data"])
      {
	GridS_bg_aer = RegularGrid<T>(Ns_bg_aer);
	GridB_bg_aer = RegularGrid<T>(Nb_bg_aer*this->Ncomposition_aer); 
	
	Background_aer_i.Resize(GridS_bg_aer, GridB_bg_aer, this->GridST2D);
	Background_aer_f.Resize(GridS_bg_aer, GridB_bg_aer, this->GridST2D);
	FileBackground_aer_i.Resize(GridS_bg_aer, GridB_bg_aer, this->GridST2D);
	FileBackground_aer_f.Resize(GridS_bg_aer, GridB_bg_aer, this->GridST2D);
	NumberBackground_aer_i.Resize(GridB_bg_aer, this->GridST2D);
	NumberBackground_aer_f.Resize(GridB_bg_aer, this->GridST2D);
	FileNumberBackground_aer_i.Resize(GridB_bg_aer, this->GridST2D);
	FileNumberBackground_aer_f.Resize(GridB_bg_aer, this->GridST2D);
    
	Background_aer_i.SetZero();
	Background_aer_f.SetZero();
	FileBackground_aer_i.SetZero();
	FileBackground_aer_f.SetZero();
	NumberBackground_aer_i.SetZero();
	NumberBackground_aer_f.SetZero();
	FileNumberBackground_aer_i.SetZero();
	FileNumberBackground_aer_f.SetZero();  
      }  

    /*** Initial concentration ***/
    if (ic_format == "Internal" &&
        this->option_process["with_external_composition"])
      {
        StreetConcentration_aer_i.Resize(GridS3D_aer, GridB3D_aer_i, GridST3D);
        StreetConcentration_aer_i.SetZero();
        StreetNumberConcentration_i.Resize(GridB2D_aer_i, this->GridST2D);
        StreetNumberConcentration_i.SetZero();
      }    
    
    /*** Aerosol street concentration ***/

    StreetConcentration_aer.Resize(GridS2D_aer, GridB2D_aer, this->GridST2D);  
    StreetNumberConcentration.Resize(GridB2D_aer, this->GridST2D);
    StreetConcentration_aer.SetZero();     
    StreetNumberConcentration.SetZero();     

    /*** Aerosol deposition ***/
    StreetSurfaceDepositedMass_aer.Resize(GridS2D_aer, GridB2D_aer, this->GridST2D);
    StreetSurfaceDepositedMass_aer.SetZero();
    StreetResuspensionRate_aer.Resize(GridS2D_aer, GridB2D_aer, this->GridST2D);
    StreetResuspensionRate_aer.SetZero();
    StreetWashoffRate_aer.Resize(GridS2D_aer, GridB2D_aer, this->GridST2D);
    StreetWashoffRate_aer.SetZero();
    StreetSurfaceDryDepositionRate_aer.Resize(GridS2D_aer, GridB2D_aer, this->GridST2D);
    StreetSurfaceDryDepositionRate_aer.SetZero();
    
    GridB_dep_aer = RegularGrid<T>(Nbin_dep_aer);

    StreetDryDepositionFlux_aer.Resize(GridB_dep_aer, this->GridST2D);
    WallDryDepositionFlux_aer.Resize(GridB_dep_aer, this->GridST2D);
    TreeDryDepositionFlux_aer.Resize(GridB_dep_aer, this->GridST2D);
    StreetDryDepositionRate_aer.Resize(GridB_dep_aer, this->GridST2D); // ug/s
    WallDryDepositionRate_aer.Resize(GridB_dep_aer, this->GridST2D);; // ug/s
    TreeDryDepositionRate_aer.Resize(GridB_dep_aer, this->GridST2D); // ug/s
       
    StreetDryDepositionFlux_aer.SetZero();
    WallDryDepositionFlux_aer.SetZero();
    TreeDryDepositionFlux_aer.SetZero();
    StreetDryDepositionRate_aer.SetZero();
    WallDryDepositionRate_aer.SetZero();
    TreeDryDepositionRate_aer.SetZero();

    /*** Scavenging ***/

    GridB_scav_aer = RegularGrid<T>(Nb_scav_aer);

    StreetScavengingFlux_aer.Resize(GridB_scav_aer, this->GridST2D);
    StreetScavengingRate_aer.Resize(GridB_scav_aer, this->GridST2D);

    StreetScavengingFlux_aer.SetZero();
    StreetScavengingRate_aer.SetZero();
    
    pH.Resize(this->GridST2D);
    // Default cloud pH.
    pH.Fill(4.5);
    SetpH();

    /*** Road traffic ***/
    RoadTraffic_2R_i.Resize(this->GridST2D);
    RoadTraffic_2R_f.Resize(this->GridST2D);
    FileRoadTraffic_2R_i.Resize(this->GridST2D);
    FileRoadTraffic_2R_f.Resize(this->GridST2D);
	
    RoadTraffic_2R_i.SetZero();
    RoadTraffic_2R_f.SetZero();
    FileRoadTraffic_2R_i.SetZero();
    FileRoadTraffic_2R_f.SetZero();

    RoadTraffic_HDV_i.Resize(this->GridST2D);
    RoadTraffic_HDV_f.Resize(this->GridST2D);
    FileRoadTraffic_HDV_i.Resize(this->GridST2D);
    FileRoadTraffic_HDV_f.Resize(this->GridST2D);
	
    RoadTraffic_HDV_i.SetZero();
    RoadTraffic_HDV_f.SetZero();
    FileRoadTraffic_HDV_i.SetZero();
    FileRoadTraffic_HDV_f.SetZero();

    RoadTraffic_PC_i.Resize(this->GridST2D);
    RoadTraffic_PC_f.Resize(this->GridST2D);
    FileRoadTraffic_PC_i.Resize(this->GridST2D);
    FileRoadTraffic_PC_f.Resize(this->GridST2D);
	
    RoadTraffic_PC_i.SetZero();
    RoadTraffic_PC_f.SetZero();
    FileRoadTraffic_PC_i.SetZero();
    FileRoadTraffic_PC_f.SetZero();

    RoadTraffic_LCV_i.Resize(this->GridST2D);
    RoadTraffic_LCV_f.Resize(this->GridST2D);
    FileRoadTraffic_LCV_i.Resize(this->GridST2D);
    FileRoadTraffic_LCV_f.Resize(this->GridST2D);
	
    RoadTraffic_LCV_i.SetZero();
    RoadTraffic_LCV_f.SetZero();
    FileRoadTraffic_LCV_i.SetZero();
    FileRoadTraffic_LCV_f.SetZero();
  }

  //! Update pH variable.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetpH()
  { 
    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        // if (this->option_process["with_pH"])
        //   {
            //! Set the meteo data
            street->SetpH(pH(st));  
          // } 
        st += 1;    
      }
  }

  //! Performs one step forward.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::Forward()
  {

    InitMeteo();

    //! Chemistry() is called within ComputeMassBalance()
    //! when non-stationary approach is used
    ComputeMassBalance();
    
    if (this->option_process["with_stationary_hypothesis"])
      if (this->option_process["with_chemistry"])
	Chemistry();
    
    this->SetStreetConcentration();

    this->SetStreetDryDepositionVelocity();
    SetStreetSurfaceDepositedMass_aer();
    SetStreetResuspensionRate_aer();
    SetStreetSurfaceDryDepositionRate_aer();
    SetStreetWashoffRate_aer();
    
    SetStreetOutput_aer();

    this->AddTime(this->Delta_t);
    this->step++;
  }

  //! Calls the functions for the transport.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetSurfaceDepositedMass_aer()
  {
    int sst = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      StreetSurfaceDepositedMass_aer(s, b, sst) =
                street->GetStreetSurfaceDepositedMass_aer(s, b);
	    }
	sst ++;
      }
  }

  //! Calls the functions for the transport.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetResuspensionRate_aer()
  {
    int sst = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      StreetResuspensionRate_aer(s, b, sst) =
                street->GetStreetSurfaceDepositedMass_aer(s, b) *
                street->GetStreetResuspensionFactor();
	    }
	sst ++;
      }
  }

  //! Calls the functions for the transport.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetSurfaceDryDepositionRate_aer()
  {
    int sst = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;
	T street_area = street->GetWidth() * street->GetLength(); // m2
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      StreetSurfaceDryDepositionRate_aer(s, b, sst) = street_area *
                street->GetStreetDryDepositionVelocity_aer(b) *
                street->GetStreetConcentration_aer(s, b);
	    }
	sst ++;
      }
  }

  //! Calls the functions for the transport.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetWashoffRate_aer()
  {
    int sst = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      StreetWashoffRate_aer(s, b, sst) =
                street->GetStreetSurfaceDepositedMass_aer(s, b) *
                street->GetStreetWashoffFactor(s);
	    }
	sst ++;
      }
  }
  
  //! Calls the functions for the transport.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::ComputeStreetSurfaceDepositedMass_aer()
  {
    int sst = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	int id_st = this->id_street(sst);
	T street_area = street->GetWidth() * street->GetLength(); // m2
	T street_volume = street->GetHeight() * street->GetWidth() *
          street->GetLength(); // m3
	T f_resusp = street->GetStreetResuspensionFactor();
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      T f_washoff = street->GetStreetWashoffFactor(s);
	      T f = f_resusp + f_washoff;
	      T street_dry_deposition_rate_aer = street_area * 
		street->GetStreetDryDepositionVelocity_aer(b); // m3/s
	      T street_scavenging_flux = street_volume *
		street->GetStreetScavengingCoefficient_aer(b);
	      T street_conc = street->GetStreetConcentration_aer(s, b);		
	      T street_surface_deposited_mass_aer =
                street->GetStreetSurfaceDepositedMass_aer(s, b);
	      if(f > 0.0)
		{
		  T D = (street_dry_deposition_rate_aer + street_scavenging_flux) *
                    street_conc; //mug/s
		  T Dmass_diff = street->GetStreetSurfaceDepositedMass_aer(s, b) - D/f;
		  street_surface_deposited_mass_aer = (D/f) +
                    Dmass_diff * exp(-f*this->Delta_t); //mug
		}
	      street->SetStreetSurfaceDepositedMass_aer(street_surface_deposited_mass_aer,
                                                        s, b);
	    }
	sst += 1;
      }
  }

  
  //! Calls the functions for the transport.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::InitMeteo()    
  {
    StreetNetworkTransport<T>::InitMeteo();
    
    if (this->option_process["with_deposition_aer"])
      ComputeDryDepositionVelocities_aer();
    if (this->option_process["with_drainage_aer"])
      CalculWashoffFactor();  
    if (this->option_process["with_scavenging_aer"])
      ComputeScavengingCoefficient_aer();
    if(this->option_process["with_resuspension"])
      CalculAerosolResuspensionFactor();
    	
  }

  //! Compute mass balance.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::ComputeMassBalance()    
  {
   
    if (this->option_process["with_stationary_hypothesis"])
      {
        if (this->option_process["with_transport"])
          {
            int niter = 0;
            const int niter_max = 1000;
            while (niter < niter_max and (!this->is_stationary))
              {
                InitInflowRate();    
                // cout << " ----> Iteration No " << niter << endl;
                ComputeInflowRateExtended();
                //! Compute the concentrations in the street-canyon.
                ComputeStreetConcentration();
                this->IsStationary(this->is_stationary);
                ++niter;

              }
            if (!this->is_stationary)
              throw string("Error: stationarity is not achieved. ") +
                "Please increase the number of iterations.";
          }
        
	if (this->option_process["with_deposition_aer"] and
            this->option_process["with_resuspension"])
	  ComputeStreetSurfaceDepositedMass_aer();	
      }
    else //no stationary
      {
	InitInflowRate();
	ComputeInflowRateExtended();
	ComputeStreetConcentrationNoStationary();  
      }
  }


  
  //! Set new street concentration and deposition
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::CalculNumberResuspension()
  {
    int b, s;
    T TotalMass, Rho_aer, Rho_tmp, NumberResuspension;
    Data<T, 1> MeanDiameter(this->Nbin_aer);
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this->Ns_aer);

    vector<T> Rho_species; // g.cm-3 -> 10-6 ug.um-3
    
    for (int s = 0; s < this->Ns_aer; s++)
      Rho_species.push_back(ssh_mass_density_layers(s));
    
    //! Compute mean diameter of each section
    MeanDiameter = 0.0;
    for(b = 0; b < this->Nbin_aer; b++)
      MeanDiameter(b) = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um
    
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	T street_volume = street->GetHeight() * street->GetWidth() *
          street->GetLength(); // m3
	for(b = 0; b < this->Nbin_aer; b++)
	  {
	    NumberResuspension = 0.0;
	    TotalMass = 0.0;
	    Rho_aer = 0.0;
	    Conc_aer_tmp.SetZero();
	    for(s = 0; s < this->Ns_aer; s++)
	      {
		TotalMass += street->GetStreetResuspension(s, b) *
                  this->Delta_t / street_volume; //mug/m3
		Conc_aer_tmp(s) = street->GetStreetResuspension(s, b) *
                  this->Delta_t / street_volume; //mug/m3
	      }
            Rho_aer = 1e9 * ComputeDensity(Conc_aer_tmp, Rho_species,
                                           TotalMass, this->Ns_aer);
 
	    NumberResuspension = TotalMass / Rho_aer / pi * 6.
	      / (MeanDiameter(b) * MeanDiameter(b) * MeanDiameter(b));
	    if (NumberResuspension * 1. != NumberResuspension)
	      throw string("NaN in NumberResuspension");
	    street->SetStreetNumberResuspension(NumberResuspension, b);
	    
	  }
      }
  }

  //Set new street concentration and deposition
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::CalculWashoffFactor()
  {
    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	T rain = street->GetRain(); //mm/h
	T groad = rain / 3600. * this->Delta_t;
	T groad_drain_min = min_water_drainage; //mm
	for (int s = 0; s < this->Ns_aer; ++s)
	  {
	    T f_washoff = 0.0;
	    T hdrain_eff = drainage_efficiency[this->species_list_aer[s]];
	    if(groad > groad_drain_min)
	      {
		f_washoff = (1 -
                             exp(-1*hdrain_eff *
                                 ((groad - groad_drain_min)/groad_drain_min))
                             ) / this->Delta_t; //s-1
		if (f_washoff * 1. != f_washoff)
		  throw string("Error ! NaN in f_washoff");
		f_washoff = max(f_washoff, 0.0);
	      }
	    street->SetStreetWashoffFactor(f_washoff, s);
	  }	
      }    
  }
  
  
  //Set new street concentration and deposition
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::CalculAerosolResuspensionFactor()
  {
    T u_ref_susp = 50.0; //km/h
    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;
	T f_resusp = 0.0;
	T rain = street->GetRain();
	if(rain < max_rain)
	  {
	    T traffic_2R = street->GetStreetRoadTraffic_2R();
	    T traffic_hdv = street->GetStreetRoadTraffic_HDV();
	    T traffic_pc = street->GetStreetRoadTraffic_PC();
	    T traffic_lcv = street->GetStreetRoadTraffic_LCV();
	    
	    T speed_2R, speed_HDV, speed_PC, speed_LCV;
	    
	    if(street->GetStreetTypo() == 1) //highway typo
	      {
		speed_2R = mean_speed_highway_2R;
		speed_HDV = mean_speed_highway_HDV;
		speed_PC = mean_speed_highway_PC;
		speed_LCV = mean_speed_highway_LCV;
	      }
	    else
	      {
		speed_2R = mean_speed_2R;
		speed_HDV = mean_speed_HDV;
		speed_PC = mean_speed_PC;
		speed_LCV = mean_speed_LCV;
	      }
	    
	    T f_resusp_2R = traffic_2R/(street->GetLength()/1000.) *
              (speed_2R/u_ref_susp) * f0_2R; //s-1
	    T f_resusp_hdv = traffic_hdv/(street->GetLength()/1000.) *
              (speed_HDV/u_ref_susp) * f0_hdv; //s-1
	    T f_resusp_pc = traffic_pc/(street->GetLength()/1000.) *
              (speed_PC/u_ref_susp) * f0_pc; //s-1
	    T f_resusp_lcv = traffic_lcv/(street->GetLength()/1000.) *
              (speed_LCV/u_ref_susp) * f0_lcv; //s-1

	    f_resusp = f_resusp_2R + f_resusp_hdv + f_resusp_pc + f_resusp_lcv; //s-1
	    
	    if (f_resusp * 1.0 != f_resusp)
	      throw string("Error! f_resusp = NaN");
	  }
	street->SetStreetResuspensionFactor(f_resusp);    
	st ++;
      }
  }
  
  //! Set new street concentration and deposition
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetOutput_aer()
  {
    SetStreetConcentration_aer();
    SetStreetNumberConcentration_aer();
    if(this->option_process["with_deposition_aer"])
      {
	SetStreetDryDepositionFlux_aer();
        SetWallDryDepositionFlux_aer();
        SetTreeDryDepositionFlux_aer();
	SetStreetDryDepositionRate_aer(); // ug/s
	SetWallDryDepositionRate_aer(); // ug/s
	SetTreeDryDepositionRate_aer(); // ug/s
      }
    if(this->option_process["with_scavenging_aer"])
      {
	SetStreetScavengingFlux_aer();
	SetStreetScavengingRate_aer();
      }
  }

  //! Sets the dry deposition flux for the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetDryDepositionFlux_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for(int b = 0; b < this->Nbin_aer; ++b)
	  for (int j = 0; j < Nbin_dep_aer; ++j)
	    if (b == j)
		StreetDryDepositionFlux_aer(j, ist) =
                  street->GetStreetDryDepositionFlux_aer(b);
        ++ist;
      }
  }

  //! Sets the wall dry deposition flux for all species at the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetWallDryDepositionFlux_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for(int b = 0; b < this->Nbin_aer; ++b)
	  for (int j = 0; j < Nbin_dep_aer; ++j)
	    if (b == j)
	      {
		WallDryDepositionFlux_aer(j, ist) =
                  street->GetWallDryDepositionFlux_aer(b);
	      }
        ++ist;
      }
  }

  //! Sets the tree leaves dry deposition flux for all species at the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetTreeDryDepositionFlux_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for(int b = 0; b < this->Nbin_aer; ++b)
	  for (int j = 0; j < Nbin_dep_aer; ++j)
	    if (b == j)
	      {
		TreeDryDepositionFlux_aer(j, ist) =
                  street->GetTreeDryDepositionFlux_aer(b);
	      }
        ++ist;
      }
  }


  //! Sets the dry deposition rate for the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetDryDepositionRate_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for(int b = 0; b < this->Nbin_aer; ++b)
	  for (int j = 0; j < Nbin_dep_aer; ++j)
	    if (b == j)
	      {
		StreetDryDepositionRate_aer(j, ist) =
                  street->GetStreetDryDepositionFlux_aer(b) *
                  street->GetLength() * street->GetWidth();
	      }
        ++ist;
      }
  }
  
  //! Sets the wall dry deposition rate for all species at the whole street-network
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetWallDryDepositionRate_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for(int b = 0; b < this->Nbin_aer; ++b)
	  for (int j = 0; j < Nbin_dep_aer; ++j)
	    if (b == j)
	      {
		WallDryDepositionRate_aer(j, ist) =
                  street->GetWallDryDepositionFlux_aer(b) *
                  street->GetHeight() * street->GetLength() * 2.0;
	      }
        ++ist;
      }
  }

  //! Sets the tree leaves dry deposition rate for all species at the whole street-network
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetTreeDryDepositionRate_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for(int b = 0; b < this->Nbin_aer; ++b)
	  for (int j = 0; j < Nbin_dep_aer; ++j)
	    if (b == j)
	      {
		TreeDryDepositionRate_aer(j, ist) =
                  street->GetTreeDryDepositionFlux_aer(b) *
                  street->GetTreeLAI() *
                  street->GetWidth() * street->GetLength();
	      }
        ++ist;
      }
  }


 //! Returns the index of a species with deposition velocities.
  /*!
    \param s species index in deposition velocities.
    \return The index of a species in the species list.
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>
  ::DepositionVelocityGlobalIndex_aer(int s) const
  {
    return this->GetSpeciesIndex_aer(DepositionVelocityName_aer(s));
  }

  //! Returns the name of a species with deposition velocities.
  /*!
    \param s species index in deposition velocities.
    \return The species name.
  */
  template<class T, class ClassChemistry>
  string StreetNetworkAerosol<T, ClassChemistry>::DepositionVelocityName_aer(int s) const
  {
    return species_list_dep_aer[s].first;
  }


  //! 
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::ComputeDryDepositionVelocities_aer()
  {
    int st = 0;
    Data<T, 1> Conc_aer_tmp(this->Ns_aer);
    Conc_aer_tmp.SetZero();
    T TotalMass;
    Array<T, 1> Rho_bin;
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc = this->Ncomposition_aer;
    vector<int> bin_dep_list;
    for (int i = 0; i < Nb_dep_aer; i++)
	bin_dep_list.push_back(i);

    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	//Particles density in each size bound
	Rho_bin.resize(Nbin_dep_aer); //unity kg/m3
	Rho_bin = 0.0;
	dist = 0;
	for (int b = 0; b < Nb_dep_aer; b++)
	  for (int ic = 0; ic < Nc; ic++)
	    {
	      TotalMass = 0.0;
	      for (int s = 0; s < this->Ns_aer; s++)
		{
		  it_begin = bin_dep_list.begin();
		  it_end = bin_dep_list.end();
		  pos = find(it_begin,it_end,b);
		  dist = distance(it_begin, pos);

		  if (dist < int(bin_dep_list.size()))
		    {
		      TotalMass = TotalMass +
                        street->GetStreetConcentration_aer(s, dist*Nc+ic);
		      Conc_aer_tmp(s) = street->GetStreetConcentration_aer(s, dist*Nc+ic);
		    }
		}
              Rho_bin(dist*Nc+ic) = 1e9 *
                ComputeDensity(Conc_aer_tmp, Rho_species_aer,
                               TotalMass, this->Ns_aer); // Density per bin in kg/m^3.
	    }

	// Street data
        T H = street->GetHeight();
        T W = street->GetWidth();
        T L = street->GetLength();
        // Compute wall and ground surface friction velocities
        T z1 = this->z0s; // altitude above the surface to compute the friction velocty (m)
        T ustar_city = street->GetStreetUstar();
        T ustar_surface = StreetNetworkTransport<T>::ComputeUstarAboveSurface(z1, H, W, this->z0s, ustar_city);

        T temperature_ = street->GetTemperature();
        T pressure_ = street->GetPressure();
	
	Array<T, 1> WetDiameter_dep_aer(Nbin_dep_aer);
	WetDiameter_dep_aer = 0.0;
	for (int b = 0; b < Nb_dep_aer; b++)
	  {
	    int pos_bin = dep_bin_list_aer[b];
	    for(int id=0; id< this->Ncomposition_aer; id++)
	      WetDiameter_dep_aer(b*this->Ncomposition_aer + id) = street->GetStreetWetDiameter_aer(pos_bin*this->Ncomposition_aer + id);
	  }

	// Deposition on tree leaves
        // The parameterization used to compute aerosol deposition on tree leaves is the same as
        // the one used for deposition on street ground and walls (zhang or giardina).
        T LAI, hmax, htrunk = 0.0;
        T ustar_htree = 0.0;
        T tree_dry_deposition_velocity = 0.0;
        T Rs_tree;
        if (this->option_process["with_tree_deposition"])
          {
            hmax = street->GetTreeHeight();
            htrunk = street->GetTrunkHeight();
            LAI = street->GetTreeLAI();
            if (LAI != 0.0 and hmax != 0.0 and htrunk != 0.0)
              {
                T sH = ComputeWangsH(H, W, hmax, htrunk, LAI, this->Cdt);
                T hm = htrunk + (hmax - htrunk) / 2; // middle height of the tree crown
                ustar_htree = ComputeWangUstarProfile(hm, H, W, this->z0s, sH, ustar_city, hmax, htrunk, LAI, this->Cdt);
              }
          }

	T Ts = temperature_;
	T Ps = pressure_;
        T DynamicViscosity = (pow(Ts,3.)*8.8848e-15) - (pow(Ts,2.)*3.2398e-11)
	  + (pow(Ts,1.)*6.2657e-08) + (pow(Ts,0.)*2.3543e-06);
	// Compute AirDensity = f(Temperature)
        T AirDensity = 1.293 * (273.15 / Ts);
        T KinematicViscosity = DynamicViscosity / AirDensity;

	for(int b = 0; b < Nbin_dep_aer; b++)
	  {
	    T wet_diameter = WetDiameter_dep_aer(b);
	    T ParticleDensity = Rho_bin(b);

            // Compute sedimentation velocity
            T L = ComputeMeanFreePath(Ts, Ps, KinematicViscosity);
            T Cu = ComputeCunninghamNumber(wet_diameter, L);
            T Vs = ComputeSedimentationVelocity(wet_diameter, Cu, KinematicViscosity, ParticleDensity);
            if (Vs == 0.0)
              throw string("Error: sedimentation velocity is zero.") + 
                "Please check particle diameter and density.\n" +
                " The diameter is " + to_str(wet_diameter) + " and the density is " +
                to_str(ParticleDensity) + ".";

            // Compute Schmidt number
            T Dm = ComputeMolecularDiffusivity(wet_diameter, Cu, Ts, DynamicViscosity);
            T Sc = ComputeSchmidtNumber(KinematicViscosity, Dm);
            
	    T Beta = 2.0; // Zhang model constant for all LUC
	    T Radius = 10.e-3; // for urban LUC
            T Rs;
            // Compute surface resistance according to Zhang et al., (2001)
            if (particles_dry_velocity_option == "zhang")
            {
              T Alpha = 1.5; // Zhang model constant for urban LUC
	      T Gamma = 0.56; // Zhang model constant for urban LUC
              T St = ComputeStokesNumber(ustar_surface, Radius, Vs, KinematicViscosity);
              T Rr = ComputeReboundCoefficient(St);
              T Eim = ComputeZhangImpactionEfficiency(St, Alpha, Beta);
              T Ein = ComputeZhangInterceptionEfficiency(wet_diameter, Radius);
              T Eb;
              if (brownian_diffusion_resistence_option == "zhang")
                Eb = ComputeZhangBrownianEfficiency(Sc, Gamma);
              else if (brownian_diffusion_resistence_option == "seigneur")
                Eb = ComputeSeigneurBrownianEfficiency(Sc);
              else
                throw string("Error in brownian_diffusion_resistence_option, choose 'zhang' or 'seigneur'");
              Rs = ComputeZhangSurfaceResistance(ustar_surface, Eim, Ein, Eb, Rr);

              if (this->option_process["with_tree_deposition"] and (LAI != 0.0 and hmax != 0.0 and htrunk != 0.0))
              {
                T St_tree = ComputeStokesNumber(ustar_htree, this->Radius_tree, Vs, KinematicViscosity);
                T Rr_tree = ComputeReboundCoefficient(St_tree);
                T Eim_tree = ComputeZhangImpactionEfficiency(St_tree, this->alpha_tree, Beta);
                T Ein_tree = ComputeZhangInterceptionEfficiency(wet_diameter, this->Radius_tree);
                T Eb_tree;
                if (brownian_diffusion_resistence_option == "zhang")
                  Eb_tree = ComputeZhangBrownianEfficiency(Sc, this->gamma_tree);
                else
                  throw string("Use the option 'zhang' for brownian_diffusion_resistance_option for deposition on tree leaves");
                Rs_tree = ComputeZhangSurfaceResistance(ustar_htree, Eim_tree, Ein_tree, Eb_tree, Rr_tree);
                tree_dry_deposition_velocity = ComputeVenkatranDepositionVelocity(Vs, Rs_tree);
              }
            }

            // Compute surface resistance according to Giardina and Buffa, (2018)
            else if (particles_dry_velocity_option == "giardina")
            {
              T tau = ComputeRelaxationTime(wet_diameter, Cu, DynamicViscosity, ParticleDensity);
              T Re_star = ComputeReynoldsNumber(ustar_surface, this->z0s, KinematicViscosity);
              T St = ComputeStokesNumber(ustar_surface, 0., Vs, KinematicViscosity); // In Giardina and Buffa, (2018) St number calculated with equation for smooth surfaces
              T Rr = ComputeReboundCoefficient(St);
              T Rii = ComputeGiardinaInertialImpactResistance(St, ustar_surface, Radius);
              T Rti = ComputeGiardinaTurbulentImpactResistance(ustar_surface, tau, KinematicViscosity);
              T Rdb; 
              if (brownian_diffusion_resistence_option == "giardina")
                Rdb = ComputeGiardinaBrownianDiffusionResistance(Sc, ustar_surface);
              else if (brownian_diffusion_resistence_option == "chamberlain")
                Rdb = ComputeChamberlainBrownianEfficiency(Sc, ustar_surface, Re_star);
              else
                throw string("Error in brownian_diffusion_resistence_option, choose 'giardina' or 'chamberlain'");
              Rs = ComputeGiardinaSurfaceResistance(Rii, Rti, Rdb, Rr);

              if (this->option_process["with_tree_deposition"] and (LAI != 0.0 and hmax != 0.0 and htrunk != 0.0))
              {
                T St_tree = ComputeStokesNumber(ustar_htree, 0., Vs, KinematicViscosity); // In Giardina and Buffa, (2018) St number calculated with equation for smooth surfaces
                T Rr_tree = ComputeReboundCoefficient(St_tree);
                T Rii_tree = ComputeGiardinaInertialImpactResistance(St_tree, ustar_htree, this->Radius_tree);
                T Rti_tree = ComputeGiardinaTurbulentImpactResistance(ustar_htree, tau, KinematicViscosity);
                T Rdb_tree;
                if (brownian_diffusion_resistence_option == "giardina")
                  Rdb_tree = ComputeGiardinaBrownianDiffusionResistance(Sc, ustar_htree);
                else
                  throw string("Use the option 'giardina' for brownian_diffusion_resistance_option for deposition on tree leaves");
                Rs_tree = ComputeGiardinaSurfaceResistance(Rii, Rti, Rdb, Rr);
                tree_dry_deposition_velocity = ComputeVenkatranDepositionVelocity(Vs, Rs_tree);
              }
            }
            else
              throw string("Error in particle dry deposition velocity option, choose 'zhang' or 'giardina'");
            
            if (Rs <= 0.0)
              throw string("Error in deposistion, Rs: ") + to_str(Rs);

            T street_dry_deposition_velocity = ComputeVenkatranDepositionVelocity(Vs, Rs);
            street->SetStreetDryDepositionVelocity_aer(street_dry_deposition_velocity, b);

	    // Sedimentation velocity does not impact deposition on walls
	    T wall_dry_deposition_velocity = 1. / Rs;
            street->SetWallDryDepositionVelocity_aer(wall_dry_deposition_velocity, b);

            street->SetTreeDryDepositionVelocity_aer(tree_dry_deposition_velocity, b);
          }
	st += 1;
      }
  }


  //! Returns the concentrations for the whole street-network.
  template<class T, class ClassChemistry>
  inline T StreetNetworkAerosol<T, ClassChemistry>
  ::GetStreetScavengingFlux_aer(int street_index, int b)
  {
    StreetNetworkTransport<T>::SetCurrentStreet(street_index);
    Street<T>* street = this->current_street;

    return street->GetStreetScavengingFlux_aer(b);
  }

  //! Sets the street scavening flux over a canopy for one specie in a specific street .
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::
  SetStreetScavengingFluxOverCanopy_aer(int street_index,
                                        T scavenging_flux_overcanopy, int b)
  {
    //SetCurrentStreet(street_index);
    StreetNetworkTransport<T>::SetCurrentStreet(street_index);
    Street<T>* street = this->current_street;
    street->SetStreetScavengingFluxOverCanopy_aer(scavenging_flux_overcanopy, b);
  }

  //! Returns scavening flux for all species in the whole street-network.
  template<class T, class ClassChemistry>
  inline Data<T, 2>& StreetNetworkAerosol<T, ClassChemistry>
  ::GetStreetScavengingFlux_aer()
  {
    return StreetScavengingFlux_aer;
  }
  
  //! Sets the scavening flux for all species in the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetScavengingFlux_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for(int b = 0; b < this->Nbin_aer; b++)
	  for (int j = 0; j < Nbin_scav_aer; ++j)
	    if (b == j)
	      StreetScavengingFlux_aer(b, ist) =
                street->GetStreetScavengingFlux_aer(b);
        ++ist;
      }
  }

  //! Returns the scavening rate for all species in the whole street-network.
  template<class T, class ClassChemistry>
  inline Data<T, 2>& StreetNetworkAerosol<T, ClassChemistry>
  ::GetStreetScavengingRate_aer()
  {
    return StreetScavengingRate_aer;
  }

  //! Sets the scavening rate for all species in the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetScavengingRate_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	for (int b = 0; b < this->Nbin_aer; b++)
	  for (int j = 0; j < Nbin_scav_aer; ++j)
	    if (b == j)
	      StreetScavengingRate_aer(b, ist) =
                street->GetStreetScavengingFlux_aer(b) *
                street->GetWidth() * street->GetLength();

        ++ist;
      }
  }

  //! Initializes scavenging coefficients for aerosols.
  /*! Computes below-cloud scavenging coefficients for aerosols.
    \param Temperature_ temperature (K).
    \param Pressure_ pressure (Pa).
    \param WetDiameter_aer_ (output) wet aerosol diameter (m).
    \param CloudHeight_ cloud basis height (m).
    \param Rain_ rain intensity (mm/h).
    \param ScavengingCoefficient_aer_ (output) the scavenging coefficients
    (s^{-1}).
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>
  ::ComputeScavengingCoefficient_aer()
  {
    int st = 0;
    int b, s, ic;
    Data<T, 1> Conc_aer_tmp(this->Ns_aer);
    Conc_aer_tmp.SetZero();
    T TotalMass;
    Array<T, 1> Rho_bin;
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc = this->Ncomposition_aer;
    Array<T, 1> coefficient(Nbin_scav_aer);
    Array<T, 1> wet_diameter_(this->Nbin_aer);
    Array<T, 1> wet_diameter_loc(Nbin_scav_aer);
    Array<T, 1> level(1);
    level = 2.;
    vector<int> bin_scav_list;
    for (int i = 0; i < Nb_scav_aer; i++)
	bin_scav_list.push_back(i);

    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
	//Particles density in each size bound
	Rho_bin.resize(Nbin_scav_aer);
	Rho_bin = 0.0;
	dist = 0;
	for (b = 0; b < Nb_scav_aer; b++)
	  for (ic = 0; ic < Nc; ic++)
	    {
	      TotalMass = 0.0;
	      for (s = 0; s < this->Ns_aer; s++)
		{
		  it_begin = bin_scav_list.begin();
		  it_end = bin_scav_list.end();
		  pos = find(it_begin,it_end,b);
		  dist = distance(it_begin, pos);

		  if (dist < int(bin_scav_list.size()))
		    {
		      TotalMass = TotalMass + StreetConcentration_aer(s,dist*Nc+ic,st);
		      Conc_aer_tmp(s) = StreetConcentration_aer(s,dist*Nc+ic,st);
		    }
		}
              Rho_bin(dist*Nc+ic) = 1e9 *
                ComputeDensity(Conc_aer_tmp,
                               Rho_species_aer,
                               TotalMass, this->Ns_aer); // Density per bin in kg / m^3.
	    }

        T temperature_ = street->GetTemperature();
        T pressure_ = street->GetPressure();
        T rain_ = street->GetRain();
        T cloudheight_ = street->GetHeight();
        T level_ = street->GetHeight();

	Array<T, 1> WetDiameter_scav_aer(this->Nbin_scav_aer);
	WetDiameter_scav_aer = 0.0;
	for (b = 0; b < Nb_scav_aer; b++)
	  {
	    int pos_bin = scav_bin_list_aer[b];
	    for(int id=0; id< this->Ncomposition_aer; id++)
	      WetDiameter_scav_aer(b*this->Ncomposition_aer + id) =
                street->GetStreetWetDiameter_aer(pos_bin*this->Ncomposition_aer + id);
	  } 
	
	_compute_scavenging_coefficient_aer(&Nbin_scav_aer,
					    &temperature_,
					    &pressure_,
					    &rain_,
					    WetDiameter_scav_aer.data(),
					    Rho_bin.data(),
					    &level_,
					    &cloudheight_,
					    coefficient.data());
	
	for (b = 0; b < Nbin_scav_aer; b++)
	  street->SetStreetScavengingCoefficient_aer(coefficient(b), b);
      }
  }

  //! Set initial_concentration array of each street.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetInitialStreetConcentration()
  {
    StreetNetworkTransport<T>::SetInitialStreetConcentration();
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      T street_conc_aer = street->GetStreetConcentration_aer(s, b); // ug/m3
	      street->SetInitialStreetConcentration_aer(street_conc_aer, s, b);
	    }
        for (int b = 0; b < this->Nbin_aer; ++b)
          {
            T street_number_conc = street->GetStreetNumberConcentration(b); // #/m3
            street->SetInitialStreetNumberConcentration(street_number_conc, b);
          }          
      }
  }
  
  
    //! (Re)Initialize the inflow rate.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::InitInflowRate()
  {
    StreetNetworkTransport<T>::InitInflowRate();
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin(); 
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
          {
            street->SetInflowRate_aer(0.0, s, b);
            street->SetMassfluxFromBackground_aer(0.0, s, b);
            street->SetMassfluxToBackground_aer(0.0, s, b);
          }
	if (this->option_process["with_number_concentration"])
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      street->SetNumberInflowRate(0.0, b);
	      street->SetNumberfluxFromBackground(0.0, b);
	      street->SetNumberfluxToBackground(0.0, b);
	    }	
      }
  }
  
  //! Compute the inflow rate with the extended matrix.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::ComputeInflowRateExtended()
  {
   
    for (typename vector<Intersection<T>* >::iterator iter = this->IntersectionVector.begin();
         iter != this->IntersectionVector.end(); iter++)
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
        Array<T, 2> extended_matrix(nstreet_inter + 1, nstreet_inter + 1);
        if (this->option_process["with_horizontal_fluctuation"])
          extended_matrix = intersection->GetGaussianMatrix();
        else
          extended_matrix = intersection->GetFluxMatrix();


        for (int j = 0; j < nstreet_inter; ++j)
          {
            int end_inter_id = 0;
            int begin_inter_id = 0;
            for (typename vector<Street<T>* >::iterator iter2 = this->StreetVector.begin();
                 iter2 != this->StreetVector.end(); iter2++)
              {
                Street<T>* street = *iter2;
                if (street_list_inter(j) == street->GetStreetID())
                  {
                    this->StreetVectorInter.push_back(street);
                    break;
                  }
              }
          }

        //! Compute outgoing volume flux (m3/s)
        for (int i = 0; i < nstreet_inter + 1; i++)
          for (int j = 0; j < nstreet_inter + 1; j++)
            if (i != j and i != 0)
              {
                T old_outgoing_flux = this->StreetVectorInter[i - 1]->GetOutgoingFlux();
                T outgoing_flux = extended_matrix(i, j) + old_outgoing_flux;
                
                this->StreetVectorInter[i - 1]->SetOutgoingFlux(outgoing_flux);
                
              }
	  

	//! Gas-phase
	for (int s = 0; s < this->Ns; ++s)
	  {
	    for (int i = 0; i < nstreet_inter + 1; i++)
	      for (int j = 0; j < nstreet_inter + 1; j++)
		{
		  if (i != j and i == 0 and j != 0)
		    {
		      T old_massflux = this->StreetVectorInter[j - 1]->GetMassfluxFromBackground(s);
		      T old_inflow_rate = this->StreetVectorInter[j - 1]->GetInflowRate(s);

		      //! Compute the inflow rate from the atmosphere to the street
		      T new_massflux  = this->StreetVectorInter[j - 1]->
                        GetBackgroundConcentration(s) * extended_matrix(i, j);
		      T massflux_from_background = new_massflux + old_massflux;
		      T inflow_rate = new_massflux + old_inflow_rate;
		      this->StreetVectorInter[j - 1]->
                        SetMassfluxFromBackground(massflux_from_background, s);
		      this->StreetVectorInter[j - 1]->
                        SetInflowRate(inflow_rate, s);
		    }
		  else if (i != j and i != 0 and j == 0)
		    {
		      T old_massflux = this->StreetVectorInter[i - 1]->
                        GetMassfluxToBackground(s);
		      //! Compute the outflow rate to the atmosphere from the street
		      T new_massflux = this->StreetVectorInter[i - 1]->
                        GetStreetConcentration(s) * extended_matrix(i, j);
		      T massflux_to_background = old_massflux + new_massflux;
		      this->StreetVectorInter[i - 1]->
                        SetMassfluxToBackground(massflux_to_background, s);
		    }
		  else if (i != j and i != 0 and j != 0)
		    {
		      T new_massflux = this->StreetVectorInter[i - 1]->
                        GetStreetConcentration(s) * extended_matrix(i, j);
		      T old_inflow_rate = this->StreetVectorInter[j - 1]->GetInflowRate(s);
		      T inflow_rate = new_massflux + old_inflow_rate;
		      this->StreetVectorInter[j - 1]->SetInflowRate(inflow_rate, s);
		    }
		}
	  }

	
	// Aerosol - mass
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      for (int i = 0; i < nstreet_inter + 1; i++)
		for (int j = 0; j < nstreet_inter + 1; j++)
		  {
		    if (i != j and i == 0 and j != 0)
		      {
			T old_massflux_aer = this->StreetVectorInter[j - 1]->
                          GetMassfluxFromBackground_aer(s, b);
			T old_inflow_rate_aer = this->StreetVectorInter[j - 1]->
                          GetInflowRate_aer(s, b);

			//! Compute the inflow rate from the atmosphere to the street
			T new_massflux_aer  = this->StreetVectorInter[j - 1]->
                          GetBackgroundConcentration_aer(s, b) * extended_matrix(i, j);
			T massflux_from_bg_aer = new_massflux_aer + old_massflux_aer;
			T inflow_rate_aer = new_massflux_aer + old_inflow_rate_aer;
			this->StreetVectorInter[j - 1]->
                          SetMassfluxFromBackground_aer(massflux_from_bg_aer, s, b);
			this->StreetVectorInter[j - 1]->
                          SetInflowRate_aer(inflow_rate_aer, s, b);
		      }
		    else if (i != j and i != 0 and j == 0)
		      {
			T old_massflux_aer = this->StreetVectorInter[i - 1]->
                          GetMassfluxToBackground_aer(s, b);
			//! Compute the outflow rate to the atmosphere from the street
			T new_massflux_aer = this->StreetVectorInter[i - 1]->
                          GetStreetConcentration_aer(s, b) * extended_matrix(i, j);
			T massflux_to_bg_aer = old_massflux_aer + new_massflux_aer;
			this->StreetVectorInter[i - 1]->
                          SetMassfluxToBackground_aer(massflux_to_bg_aer, s, b);
		      }
		    else if (i != j and i != 0 and j != 0)
		      {
			T new_massflux_aer = this->StreetVectorInter[i - 1]->
                          GetStreetConcentration_aer(s, b) * extended_matrix(i, j);
			T old_inflow_rate_aer = this->StreetVectorInter[j - 1]->
                          GetInflowRate_aer(s, b);
			T inflow_rate_aer = new_massflux_aer + old_inflow_rate_aer;
			this->StreetVectorInter[j - 1]->
                          SetInflowRate_aer(inflow_rate_aer, s, b);
		      }
		  }
	    } 
	if (this->option_process["with_number_concentration"])
	  {
	    // Aerosol - number
	    for (int b = 0; b < this->Nbin_aer; ++b)
	      for (int i = 0; i < nstreet_inter + 1; i++)
		for (int j = 0; j < nstreet_inter + 1; j++)
		  {
		    if (i != j and i == 0 and j != 0)
		      {
			T old_numberflux = this->StreetVectorInter[j - 1]->
                          GetNumberfluxFromBackground(b);
			T old_number_inflow_rate = this->StreetVectorInter[j - 1]->
                          GetNumberInflowRate(b);

			//! Compute the inflow rate from the atmosphere to the street
			T new_numberflux  = this->StreetVectorInter[j - 1]->
                          GetBackgroundNumberConcentration(b) * extended_matrix(i, j);
			T numberflux_from_bg = new_numberflux + old_numberflux;
			T number_inflow_rate = new_numberflux + old_number_inflow_rate;
			this->StreetVectorInter[j - 1]->
                          SetNumberfluxFromBackground(numberflux_from_bg, b);
			this->StreetVectorInter[j - 1]->
                          SetNumberInflowRate(number_inflow_rate, b);
		      }
		    else if (i != j and i != 0 and j == 0)
		      {
			T old_numberflux = this->StreetVectorInter[i - 1]->
                          GetNumberfluxToBackground(b);
			//! Compute the outflow rate to the atmosphere from the street
			T new_numberflux = this->StreetVectorInter[i - 1]->
                          GetStreetNumberConcentration(b) * extended_matrix(i, j);
			T numberflux_to_bg = old_numberflux + new_numberflux;
			this->StreetVectorInter[i - 1]->
                          SetNumberfluxToBackground(numberflux_to_bg, b);
		      }
		    else if (i != j and i != 0 and j != 0)
		      {
			T new_numberflux = this->StreetVectorInter[i - 1]->
                          GetStreetNumberConcentration(b) * extended_matrix(i, j);
			T old_number_inflow_rate = this->
                          StreetVectorInter[j - 1]->GetNumberInflowRate(b);
			T number_inflow_rate = new_numberflux + old_number_inflow_rate;
			this->StreetVectorInter[j - 1]->
                          SetNumberInflowRate(number_inflow_rate, b);
		      }
		  }
	  }

	//! Clear the street vector for the intersection. 
	this->StreetVectorInter.clear();
      }
  } 
  
  
  //! Compute the concentrations in the street-canyon using the flux balance equation.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::ComputeStreetConcentration()
  {
    const T delta_concentration_min = 0.01;
    int st = 0;

    vector<T> Rho_species; // g.cm-3 -> 10-6 ug.um-3
    if(number_computation_option == "based_on_mass")
      {
        for (int s = 0; s < this->Ns_aer; s++)
          Rho_species.push_back(ssh_mass_density_layers(s));
      }

    StreetNetworkTransport<T>::ComputeStreetConcentration();
    
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
	
        Street<T>* street = *iter;
        T transfer_velocity = street->GetTransferVelocity(); // m/s
        T temp = transfer_velocity * street->GetWidth() * street->GetLength();// m3/s
        T outgoing_flux = street->GetOutgoingFlux(); // m3/s
        T street_volume = street->GetHeight() * street->GetWidth() * street->GetLength(); // m3

        T street_area = street->GetWidth() * street->GetLength(); // m2
	
        //! symmetric walls
        T wall_area = 2.0 * street->GetHeight() * street->GetLength(); // m2

        T rain = street->GetRain();
	bool is_stationary_local = true;

   	//! Aerosol - mass
        for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    {
	      T street_conc_aer = street->GetStreetConcentration_aer(s, b); // ug/m3
	      T init_street_conc_aer = street->GetInitialStreetConcentration_aer(s, b); // ug/m3
	      
	      T emission_rate_aer = street->GetEmission_aer(s, b); //  ug/s
	      T inflow_rate_aer = street->GetInflowRate_aer(s, b); // ug/s
	      T conc_bg_aer = street->GetBackgroundConcentration_aer(s, b); // ug/m3

	      T deposition_flow_rate_aer = 0.0;
	      T scavenging_flow_rate_aer = 0.0;
	      T resuspension_rate = 0.0;
              
	      if(this->option_process["with_deposition_aer"])
		{
		  T street_dry_deposition_flow_rate_aer = street_area * 
		    street->GetStreetDryDepositionVelocity_aer(b); // m3/s
		  T wall_dry_deposition_flow_rate_aer = wall_area * 
		    street->GetWallDryDepositionVelocity_aer(b); // m3/s

                  if (this->option_process["with_tree_deposition"])
                    {
                      T leaf_surface = street->GetTreeLAI() * street_area; // m2
                      T tree_dry_deposition_flow_rate_aer = leaf_surface *
                        street->GetTreeDryDepositionVelocity_aer(b); // m3/s
                      deposition_flow_rate_aer = street_dry_deposition_flow_rate_aer +
                        wall_dry_deposition_flow_rate_aer + tree_dry_deposition_flow_rate_aer; // m3/s
                    }
                  else
                    deposition_flow_rate_aer = street_dry_deposition_flow_rate_aer +
                      wall_dry_deposition_flow_rate_aer; // m3/s

                  
		  if(this->option_process["with_resuspension"])
		    resuspension_rate = street->GetStreetResuspensionFactor() *
                      street->GetStreetSurfaceDepositedMass_aer(s, b);		    
		}
	      if(this->option_process["with_scavenging_aer"])
		if(rain > 0.0)
		  scavenging_flow_rate_aer = street_volume * 
		    street->GetStreetScavengingCoefficient_aer(b); // m3/s
	      
	      //! Compute the new concentrations.
	      T street_conc_new_aer = 0.0;
	      if ((temp + outgoing_flux + deposition_flow_rate_aer + scavenging_flow_rate_aer) == 0.0)
		{
		  //! The total outgoing flux should not be zero
		  //! because it leads to an infinite increase of pollutant concentrations
		  //! in Munich which uses a steady-state approximation.
		  //! In this case, the concentrations are assumed to be constant
		  //! because the total incoming flux should be equal to zero.
		  street_conc_new_aer = ((emission_rate_aer + resuspension_rate) * this->Delta_t) / street_volume
		    + init_street_conc_aer;
		}
	      else
		{
		  street_conc_new_aer = (emission_rate_aer + inflow_rate_aer + temp * conc_bg_aer + resuspension_rate) /
		    (temp + outgoing_flux + deposition_flow_rate_aer + scavenging_flow_rate_aer);
		}

	      //! Set the minimum possible concentrations.
	      street_conc_new_aer = max(street_conc_new_aer, 0.0);
	      
	      //! Check the stationarity
	      T delta_concentration_aer = abs(street_conc_new_aer - street_conc_aer);
	      if ((street_conc_aer != 0.0) and
                  (delta_concentration_aer > delta_concentration_min))
		is_stationary_local = false;

	      //! Set the new concentrations.
	      street->SetStreetConcentration_aer(street_conc_new_aer, s, b);

	      T massflux_roof_to_bg_aer;
	      if (is_stationary_local)
		{
		  massflux_roof_to_bg_aer = temp *
                    (street_conc_new_aer - conc_bg_aer); // ug/s
		  street->SetMassfluxRoofToBackground_aer(massflux_roof_to_bg_aer, s, b);		  

		  T conc_delta_aer = street_conc_new_aer - init_street_conc_aer;
		  T street_quantity_delta_aer = conc_delta_aer * street_volume; // ug
		  street->SetStreetQuantityDelta_aer(street_quantity_delta_aer, s, b);
		}
	      street->SetStreetDryDepositionFlux_aer
                (street_conc_new_aer *
                 street->GetStreetDryDepositionVelocity_aer(b),b); // ug/m2/s
            
	      street->SetWallDryDepositionFlux_aer
                (street_conc_new_aer * 
                 street->GetWallDryDepositionVelocity_aer(b), b); // ug/m2/s

	      street->SetTreeDryDepositionFlux_aer
                (street_conc_new_aer *
                 street->GetTreeDryDepositionVelocity_aer(b),b); // ug/m2/s

	      T street_scavenging_flux_incanopy_aer = street_conc_new_aer * 
		street->GetStreetScavengingCoefficient_aer(b) *
		street->GetHeight();

	      T street_scavenging_flux_total_aer = street_scavenging_flux_incanopy_aer +
                street->GetStreetScavengingFluxOverCanopy_aer(b);
	      street->SetStreetScavengingFlux_aer(street_scavenging_flux_total_aer, b); // ug/m2/s
	    
	    }

        if (this->option_process["with_number_concentration"])
	  {
	    for (int b = 0; b < this->Nbin_aer; ++b)
	      {
		T street_number_conc = street->GetStreetNumberConcentration(b); // #/m3
		T street_number_conc_new = 0.0;
		T delta_concentration_min_aer = 0.01*street_number_conc;
		T number_conc_bg = street->GetBackgroundNumberConcentration(b); // #/m3
		T init_street_number_conc = street->GetInitialStreetNumberConcentration(b); // #/m3
		    
		if (number_computation_option == "based_on_transport")
		  {
		    //Aerosol - number
				  
		    T number_emission_rate = street->GetNumberEmission(b); //  #/s
		    T number_inflow_rate = street->GetNumberInflowRate(b); // #/s
		    T number_deposition_flow_rate = 0.0;
		    T number_scavenging_flow_rate = 0.0;
		    T number_resuspension_rate = 0.0;

		    if(this->option_process["with_deposition_aer"])
		      {
			T street_number_dry_deposition_flow_rate = street_area * 
			  street->GetStreetDryDepositionVelocity_aer(b); // m3/s
			T wall_number_dry_deposition_flow_rate = wall_area * 
			  street->GetWallDryDepositionVelocity_aer(b); // m3/s

                        if (this->option_process["with_tree_deposition"])
                          {
                            T tree_number_dry_deposition_flow_rate = street->GetTreeLAI() * street_area * 
                              street->GetTreeDryDepositionVelocity_aer(b); // m3/s
                            number_deposition_flow_rate = street_number_dry_deposition_flow_rate +
                              wall_number_dry_deposition_flow_rate + tree_number_dry_deposition_flow_rate; // m3/s
                          }
                        else
                          number_deposition_flow_rate = street_number_dry_deposition_flow_rate +
                            wall_number_dry_deposition_flow_rate; // m3/s

			if(this->option_process["with_resuspension"])
			  number_resuspension_rate = street->GetStreetResuspensionFactor() * street->GetStreetSurfaceDepositedNumber(b);
		      }
	      
		    if(this->option_process["with_scavenging_aer"])
		      if(rain > 0.0)
			number_scavenging_flow_rate = street_volume * 
			  street->GetStreetScavengingCoefficient(b); // m3/s
		

		    //! Compute the new concentrations.
		    street_number_conc_new = 0.0;
		    if ((temp + outgoing_flux + number_deposition_flow_rate + number_scavenging_flow_rate) == 0.0)
		      {
			//! The total outgoing flux should not be zero
			//! because it leads to an infinite increase of pollutant concentrations
			//! in Munich which uses a steady-state approximation.
			//! In this case, the concentrations are assumed to be constant
			//! because the total incoming flux should be equal to zero.

			street_number_conc_new = (number_emission_rate * this->Delta_t) / street_volume + init_street_number_conc;
		      }
		    else
		      street_number_conc_new = (number_emission_rate + number_inflow_rate + temp * number_conc_bg) /
			(temp + outgoing_flux + number_deposition_flow_rate + number_scavenging_flow_rate);
		  }

		if (number_computation_option == "based_on_mass")
		  {
		    Data<T, 1> concentration_aer_bin(this->Ns_aer);
		    concentration_aer_bin.SetZero();
		    T TotalMass, Rho_aer;
		    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um
		    TotalMass = 0.0;
		    for (int s = 0; s < this->Ns_aer - 1; s++)
		      {
			concentration_aer_bin(s) = street->GetStreetConcentration_aer(s, b);
			TotalMass += concentration_aer_bin(s);
		      }

                    if (this->option_process["with_fixed_density"])
                      // conversion from kg/m3 to µg/µm3
                      Rho_aer = fixed_density_aer * 1.e-9;
                    else
                      Rho_aer = ComputeDensity(concentration_aer_bin, Rho_species, TotalMass, this->Ns_aer);
		    street_number_conc_new = TotalMass/Rho_aer/pi*6./
                      (MeanDiameter*MeanDiameter*MeanDiameter);
		  }

		//! Set the minimum possible concentrations.
		street_number_conc_new = max(street_number_conc_new, 0.0);

		//! Check the stationarity
		T delta_number_concentration = abs(street_number_conc_new - street_number_conc);

		if ((street_number_conc != 0.0) and
                    (delta_number_concentration > delta_concentration_min_aer))
		  is_stationary_local = false;
		   
		//! Set the new concentrations.
		street->SetStreetNumberConcentration(street_number_conc_new, b);

		T numberflux_roof_to_bg;
		if (is_stationary_local)
		  {
		    numberflux_roof_to_bg = temp * (street_number_conc_new - number_conc_bg); // #/s
		    street->SetNumberfluxRoofToBackground(numberflux_roof_to_bg, b);

		    T number_conc_delta = street_number_conc_new - init_street_number_conc;
		    T street_number_quantity_delta = number_conc_delta * street_volume; // #
		    street->SetStreetNumberQuantityDelta(street_number_quantity_delta, b);
		  }

		street->SetStreetNumberDryDepositionFlux(street_number_conc_new * 
							 street->GetStreetDryDepositionVelocity_aer(b), b); // ug/m2/s
            
		street->SetWallNumberDryDepositionFlux(street_number_conc_new * 
						       street->GetWallDryDepositionVelocity_aer(b), b); // ug/m2/s

		street->SetTreeNumberDryDepositionFlux(street_number_conc_new * 
							 street->GetTreeDryDepositionVelocity_aer(b), b); // ug/m2/s

		T street_number_scavenging_flux_incanopy = street_number_conc_new * 
		  street->GetStreetScavengingCoefficient_aer(b) *
		  street->GetHeight();

		T street_number_scavenging_flux_total = street_number_scavenging_flux_incanopy + street->GetStreetNumberScavengingFluxOverCanopy(b);
		street->SetStreetNumberScavengingFlux(street_number_scavenging_flux_total, b); // ug/m2/s

		  
	      }
	  }

        //! True if stationarity is obtained for a street.
        street->SetStationary(is_stationary_local);
	st += 1;
      } 
  }

  //! Compute the concentrations in the street-canyon using the flux balance equation.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::ComputeStreetConcentrationNoStationary()
  {
    int zero = 0;
    T sub_delta_t_init, sub_delta_t;
    sub_delta_t_init = 0.0;
    sub_delta_t = 0.0;
    
    vector<T> Rho_species; // g.cm-3 -> 10-6 ug.um-3
    if(number_computation_option == "based_on_mass")
      {
        for (int s = 0; s < this->Ns_aer; s++)
          Rho_species.push_back(ssh_mass_density_layers(s));
      }

    // MPI implementation.
    // 'For' loop on the street segments is parallelized.
    int first_index_along_source, last_index_along_source;
    
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    
    int nstreet = this->StreetVector.size();

    // The following arrays need to be used only for MPI.
    Array<T, 2> concentration_mpi;
    Array<T, 2> massflux_roof_to_background_mpi;
    Array<T, 2> street_quantity_delta_mpi;

    Array<T, 3> concentration_aer_mpi;
    Array<T, 3> massflux_roof_to_background_aer_mpi;
    Array<T, 3> street_quantity_delta_aer_mpi;
    Array<T, 2> concentration_number_mpi;
    
    concentration_mpi.resize(nstreet, this->Ns);
    massflux_roof_to_background_mpi.resize(nstreet, this->Ns);
    street_quantity_delta_mpi.resize(nstreet, this->Ns);

    concentration_aer_mpi.resize(nstreet, this->Ns_aer, this->Nbin_aer);
    massflux_roof_to_background_aer_mpi.resize(nstreet,
                                               this->Ns_aer, this->Nbin_aer);
    street_quantity_delta_aer_mpi.resize(nstreet,
                                         this->Ns_aer, this->Nbin_aer);
    concentration_number_mpi.resize(nstreet, this->Nbin_aer);
    
    concentration_mpi = 0.0;
    massflux_roof_to_background_mpi = 0.0;
    street_quantity_delta_mpi = 0.0;

    concentration_aer_mpi= 0.0;
    massflux_roof_to_background_aer_mpi = 0.0;
    street_quantity_delta_aer_mpi = 0.0;
    concentration_number_mpi = 0.0;
    
    // Scatter buffers for Street-type arrays to all processes.
    BaseModuleParallel::ScatterSlice_source_MPI(concentration_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(massflux_roof_to_background_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(street_quantity_delta_mpi);

    BaseModuleParallel::ScatterSlice_source_MPI(concentration_aer_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(massflux_roof_to_background_aer_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(street_quantity_delta_aer_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(concentration_number_mpi);
    
    // Get the first and last indices of streets for each process.
    // The number of slices is the number of streets / the number of processes.
    // For example, if 1000 streets are computed by 10 processes,
    // Process 1: Street 1 (first_index) to 100 (last_index),
    // Process 2: Street 101 (first_index) to 200 (last_index), ...
    // See BaseModuleParallel.cxx
    BaseModuleParallel::GetEdgePartition_source
      (first_index_along_source, last_index_along_source);

#else
    first_index_along_source = 0;
    last_index_along_source = this->StreetVector.size();
#endif
  
    // Parallelized loop.
    for (int i = first_index_along_source; i < last_index_along_source; i++)
      {

        Street<T>* street = this->StreetVector.at(i);

        T transfer_velocity = street->GetTransferVelocity(); // m/s
        T temp = transfer_velocity * street->GetWidth() * street->GetLength();// m3/s
        T outgoing_flux = street->GetOutgoingFlux(); // m3/s
        T inflow_flux = street->GetIncomingFlux();
	T street_volume = street->GetHeight() * street->GetWidth() * street->GetLength(); // m3
        T street_area = street->GetWidth() * street->GetLength(); // m2
        //! symmetric walls
        T wall_area = 2.0 * street->GetHeight() * street->GetLength(); // m2
        T rain = street->GetRain();
	
	//gas phase
	Array<T, 1> concentration_array(this->Ns);
	Array<T, 1> concentration_array_tmp(this->Ns);
	Array<T, 1> init_concentration_array(this->Ns);
	Array<T, 1> background_concentration_array(this->Ns);
	Array<T, 1> new_concentration_array(this->Ns);
	Array<T, 1> emission_rate_array(this->Ns);
	Array<T, 1> inflow_rate_array(this->Ns);
	Array<T, 1> deposition_flow_rate_array(this->Ns);
	Array<T, 1> street_deposition_flow_rate_array(this->Ns);
	Array<T, 1> scavenging_flow_rate_array(this->Ns);
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
	deposition_flow_rate_array = 0.0;
	street_deposition_flow_rate_array = 0.0;
	scavenging_flow_rate_array = 0.0;
	street_surface_deposited_mass_array = 0.0;
	new_street_surface_deposited_mass_array = 0.0;
	washoff_factor_array = 0.0;
	
	for (int s = 0; s < this->Ns; ++s)
          {
	    concentration_array(s) = street->GetStreetConcentration(s);
	    init_concentration_array(s) = street->GetStreetConcentration(s);
	    background_concentration_array(s) = street->GetBackgroundConcentration(s);
	    emission_rate_array(s) = street->GetEmission(s);
	    inflow_rate_array(s) = street->GetInflowRate(s);

	    if(this->option_process["with_deposition"])
	      {
		street_deposition_flow_rate_array(s) = street_area * 
		  street->GetStreetDryDepositionVelocity(s); // m3/s
		T wall_dry_deposition_flow_rate = wall_area * 
		  street->GetWallDryDepositionVelocity(s); // m3/s

                if (this->option_process["with_tree_deposition"])
                  {
                    T tree_radius = (street->GetTreeHeight() - street->GetTrunkHeight())/2.; // m
                    T leaf_surface = street->GetTreeLAI() * street_area; // m2
                    T tree_dry_deposition_flow_rate = leaf_surface *
                      street->GetTreeDryDepositionVelocity(s); // m3/s
                    deposition_flow_rate_array(s) = street_deposition_flow_rate_array(s) +
                      wall_dry_deposition_flow_rate + tree_dry_deposition_flow_rate; // m3/s
                  }
                else
                  deposition_flow_rate_array(s) = street_deposition_flow_rate_array(s) +
                      wall_dry_deposition_flow_rate; // m3/s
                  
	      }
	    
	    if(this->option_process["with_scavenging"])
	      if(rain > 0.0)
		scavenging_flow_rate_array(s) = street_volume * 
		  street->GetStreetScavengingCoefficient(s); // m3/s
	  }

	StreetNetworkTransport<T>::
	InitStep(sub_delta_t_init,
		 this->sub_delta_t_min,
		 transfer_velocity,
		 temp,
		 outgoing_flux,
		 street_volume,
		 concentration_array,
		 background_concentration_array,
		 emission_rate_array,
		 inflow_rate_array,
		 deposition_flow_rate_array,
		 street_deposition_flow_rate_array,
		 scavenging_flow_rate_array,
		 zero, //resuspension factor
		 washoff_factor_array, //zero for gas phase
		 street_surface_deposited_mass_array, //zero for gas phase
		 this->Ns,
		 zero); //bin

	
	//aerosol mass
	T resuspension_factor = street->GetStreetResuspensionFactor();
	
	Array<T, 2> new_concentration_array_aer(this->Ns_aer, this->Nbin_aer);
	Array<T, 2> concentration_array_aer_tmp(this->Ns_aer, this->Nbin_aer);
	Array<T, 2> init_concentration_array_aer(this->Ns_aer, this->Nbin_aer);
	Array<T, 1> washoff_factor_array_aer(this->Ns_aer);
	
	new_concentration_array_aer = 0.0;
	concentration_array_aer_tmp = 0.0;
	init_concentration_array_aer = 0.0;
	washoff_factor_array_aer = 0.0;

	T sub_delta_t_init_aer = 0.0;
	Array<T, 1> sub_delta_t_init_bin(this->Nbin_aer);
	Array<T, 1> deposition_flow_rate_array_aer(this->Nbin_aer);
	Array<T, 1> street_deposition_flow_rate_array_aer(this->Nbin_aer);
	Array<T, 1> scavenging_flow_rate_array_aer(this->Nbin_aer);

	sub_delta_t_init_bin = 0.0;
	deposition_flow_rate_array_aer = 0.0;
	street_deposition_flow_rate_array_aer = 0.0;
	scavenging_flow_rate_array_aer = 0.0;

	//concentration for each species at each size bound
	Array<T, 1> init_concentration_array_bin(this->Ns_aer);
	Array<T, 1> concentration_array_bin(this->Ns_aer);
	Array<T, 1> concentration_array_bin_tmp(this->Ns_aer);
	Array<T, 1> background_concentration_array_bin(this->Ns_aer);
	Array<T, 1> emission_rate_array_bin(this->Ns_aer);
	Array<T, 1> inflow_rate_array_bin(this->Ns_aer);
	Array<T, 1> new_concentration_array_bin(this->Ns_aer);
	Array<T, 1> street_surface_deposited_mass_array_bin(this->Ns_aer);
	Array<T, 1> new_street_surface_deposited_mass_array_bin(this->Ns_aer);

	init_concentration_array_bin = 0.0;
	concentration_array_bin = 0.0;
	concentration_array_bin_tmp = 0.0;
	background_concentration_array_bin = 0.0;
	emission_rate_array_bin = 0.0;
	inflow_rate_array_bin = 0.0;
	deposition_flow_rate_array_aer = 0.0;
	street_deposition_flow_rate_array_aer = 0.0;
	scavenging_flow_rate_array_aer = 0.0;
	new_concentration_array_bin = 0.0;
	street_surface_deposited_mass_array_bin = 0.0;
	new_street_surface_deposited_mass_array_bin = 0.0;

	//number declarations
	Array<T, 1> number_concentration_array(this->Nbin_aer);
	Array<T, 1> new_number_concentration_array(this->Nbin_aer);
	Array<T, 1> number_concentration_array_tmp(this->Nbin_aer);
	Array<T, 1> background_number_concentration_array(this->Nbin_aer);
	Array<T, 1> number_emission_rate_array(this->Nbin_aer);
	Array<T, 1> number_inflow_rate_array(this->Nbin_aer);
	Array<T, 1> street_surface_deposited_number_array(this->Nbin_aer);
	Array<T, 1> new_street_surface_deposited_number_array(this->Nbin_aer);

	number_concentration_array = 0.0;
	new_number_concentration_array = 0.0;
	number_concentration_array_tmp = 0.0;
	background_number_concentration_array = 0.0;
	number_emission_rate_array = 0.0;
	number_inflow_rate_array = 0.0;
	street_surface_deposited_number_array = 0.0;
	new_street_surface_deposited_number_array = 0.0;
	
	sub_delta_t_init_bin = 0.0;
	sub_delta_t_init_aer = 0.0;
	
	
	for (int b = 0; b < this->Nbin_aer; ++b)
	  {
	    if(this->option_process["with_deposition_aer"])
	      {
		street_deposition_flow_rate_array_aer(b) = street_area * 
		  street->GetStreetDryDepositionVelocity_aer(b); // m3/s

		T wall_dry_deposition_flow_rate_aer = wall_area * 
		  street->GetWallDryDepositionVelocity_aer(b); // m3/s

                if (this->option_process["with_tree_deposition"])
                  {
                    T tree_radius = (street->GetTreeHeight() - street->GetTrunkHeight())/2.; // m
                    T leaf_surface = street->GetTreeLAI() * street_area; // m2
                    T tree_dry_deposition_flow_rate_aer = leaf_surface *
                      street->GetTreeDryDepositionVelocity_aer(b); // m3/s

                    deposition_flow_rate_array_aer(b) = street_deposition_flow_rate_array_aer(b) +
                      wall_dry_deposition_flow_rate_aer + tree_dry_deposition_flow_rate_aer; // m3/s
                  }
                else
                  deposition_flow_rate_array_aer(b) = street_deposition_flow_rate_array_aer(b) +
                    wall_dry_deposition_flow_rate_aer; // m3/s                  
	      }
	    if(this->option_process["with_scavenging_aer"])
	      if(rain > 0.0)
		scavenging_flow_rate_array_aer(b) = street_volume * 
		  street->GetStreetScavengingCoefficient_aer(b); // m3/s;
	    for(int s = 0; s < this->Ns_aer; ++s)
	      {
		init_concentration_array_aer(s, b) = street->GetStreetConcentration_aer(s,b);
		init_concentration_array_bin(s) = street->GetStreetConcentration_aer(s,b);
		concentration_array_bin(s) = street->GetStreetConcentration_aer(s,b);
		background_concentration_array_bin(s) = street->GetBackgroundConcentration_aer(s,b);
		emission_rate_array_bin(s) = street->GetEmission_aer(s,b);
		inflow_rate_array_bin(s) = street->GetInflowRate_aer(s,b);
		street_surface_deposited_mass_array_bin(s) = street->GetStreetSurfaceDepositedMass_aer(s, b);
		washoff_factor_array_aer(s) = street->GetStreetWashoffFactor(s);
	      }


            StreetNetworkTransport<T>::
	      InitStep(sub_delta_t_init_aer,
		       this->sub_delta_t_min,
		       transfer_velocity,
		       temp,
		       outgoing_flux,
		       street_volume,
		       concentration_array_bin,
		       background_concentration_array_bin,
		       emission_rate_array_bin,
		       inflow_rate_array_bin,
		       deposition_flow_rate_array_aer,
		       street_deposition_flow_rate_array_aer,
		       scavenging_flow_rate_array_aer,
		       resuspension_factor,
		       washoff_factor_array_aer,
		       street_surface_deposited_mass_array_bin,
		       this->Ns_aer,
		       b);
	    sub_delta_t_init_bin(b) = sub_delta_t_init_aer;
	  }

	sub_delta_t_init_aer = min(sub_delta_t_init_bin);
	sub_delta_t_init = min(sub_delta_t_init, sub_delta_t_init_aer);
	
	Date current_date_tmp = this->current_date;
	Date next_date = this->current_date;
	next_date.AddSeconds(this->Delta_t);
	Date next_date_tmp = this->current_date;
	next_date_tmp.AddSeconds(sub_delta_t_init);

	while (current_date_tmp < next_date)
	  {
            if (this->option_process["with_transport"])
              {
                
                //gas phase
                //! Get street concentrations.
                for (int s = 0; s < this->Ns; ++s)
                  concentration_array(s) = street->GetStreetConcentration(s);

                //! Use the ETR method to calculates new street concentrations.
                if (this->option_method == "ETR")
                  {
                    StreetNetworkTransport<T>::
                      ETRConcentration(transfer_velocity,
                                       temp,
                                       outgoing_flux,
                                       street_volume,
                                       concentration_array,
                                       concentration_array_tmp,
                                       background_concentration_array,
                                       emission_rate_array,
                                       inflow_rate_array,
                                       deposition_flow_rate_array,
                                       street_deposition_flow_rate_array,
                                       scavenging_flow_rate_array,
                                       zero, //resuspension factor
                                       washoff_factor_array, //zero for gas-phase
                                       street_surface_deposited_mass_array,
                                       new_concentration_array,
                                       sub_delta_t_init,
                                       this->Ns,
                                       zero,
                                       new_street_surface_deposited_mass_array);
                  }
                else if(this->option_method == "Rosenbrock")
                  {
                    StreetNetworkTransport<T>::
                      RosenbrockConcentration(transfer_velocity,
                                              temp,
                                              outgoing_flux,
                                              street_volume,
                                              concentration_array,
                                              concentration_array_tmp,
                                              background_concentration_array,
                                              emission_rate_array,
                                              inflow_rate_array,
                                              deposition_flow_rate_array,
                                              new_concentration_array,
                                              sub_delta_t_init,
                                              inflow_flux,
                                              street_deposition_flow_rate_array,
                                              scavenging_flow_rate_array,
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
	    
                //aerosol mass
                //Transport is calculated separately for each size bound
                for (int b = 0; b < this->Nbin_aer; ++b)
                  {
                    Array<T, 1> deposition_mass_bin(this->Ns_aer);
                    deposition_mass_bin = 0.0;
		
                    for (int s = 0; s < this->Ns_aer; ++s)
                      {
                        concentration_array_bin(s) = street->GetStreetConcentration_aer(s,b);
                        background_concentration_array_bin(s) = street->GetBackgroundConcentration_aer(s,b);
                        emission_rate_array_bin(s) = street->GetEmission_aer(s,b);
                        inflow_rate_array_bin(s) = street->GetInflowRate_aer(s,b);
                        street_surface_deposited_mass_array_bin(s) = street->GetStreetSurfaceDepositedMass_aer(s, b);
                        washoff_factor_array_aer(s) = street->GetStreetWashoffFactor(s);
                      }

                    if (this->option_method == "ETR")
                      {
                        StreetNetworkTransport<T>::
                          ETRConcentration(transfer_velocity,
                                           temp,
                                           outgoing_flux,
                                           street_volume,
                                           concentration_array_bin,
                                           concentration_array_bin_tmp,
                                           background_concentration_array_bin,
                                           emission_rate_array_bin,
                                           inflow_rate_array_bin,
                                           deposition_flow_rate_array_aer,
                                           street_deposition_flow_rate_array_aer,
                                           scavenging_flow_rate_array_aer,
                                           resuspension_factor,
                                           washoff_factor_array_aer,
                                           street_surface_deposited_mass_array_bin,
                                           new_concentration_array_bin,
                                           sub_delta_t_init,
                                           this->Ns_aer,
                                           b,
                                           new_street_surface_deposited_mass_array_bin);
                      }
                    else if(this->option_method == "Rosenbrock")
                      {
                        StreetNetworkTransport<T>::
                          RosenbrockConcentration(transfer_velocity,
                                                  temp,
                                                  outgoing_flux,
                                                  street_volume,
                                                  concentration_array_bin,
                                                  concentration_array_bin_tmp,
                                                  background_concentration_array_bin,
                                                  emission_rate_array_bin,
                                                  inflow_rate_array_bin,
                                                  deposition_flow_rate_array_aer,
                                                  new_concentration_array_bin,
                                                  sub_delta_t_init,
                                                  inflow_flux,
                                                  street_deposition_flow_rate_array_aer,
                                                  scavenging_flow_rate_array_aer,
                                                  resuspension_factor,
                                                  washoff_factor_array_aer,
                                                  street_surface_deposited_mass_array_bin,
                                                  new_street_surface_deposited_mass_array_bin,
                                                  this->Ns_aer,
                                                  b);
                      }
                    for (int s = 0; s < this->Ns_aer; ++s)
                      {
                        street->SetStreetSurfaceDepositedMass_aer(new_street_surface_deposited_mass_array_bin(s), s, b);
                        street->SetStreetConcentration_aer(new_concentration_array_bin(s), s, b);
                        concentration_array_aer_tmp(s, b) = concentration_array_bin_tmp(s);
                        new_concentration_array_aer(s, b) = new_concentration_array_bin(s);
                      }
                  }
            
                //aerosol number
                if (this->option_process["with_number_concentration"])
                  {

                    // Update the number concentration from SSH-aerosol
                    for (int b = 0; b < this->Nbin_aer; ++b)
                      {
                        number_concentration_array(b) =
                          street->GetStreetNumberConcentration(b);
                      }
                    
                    if(number_computation_option == "based_on_transport")
                      {
                        for (int b = 0; b < this->Nbin_aer; ++b)
                          {
                            street_surface_deposited_number_array(b) = street->GetStreetSurfaceDepositedNumber(b);
                            background_number_concentration_array(b) = street->GetBackgroundNumberConcentration(b); // #/m3
                            number_emission_rate_array(b) = street->GetNumberEmission(b); //  #/s
                            number_inflow_rate_array(b) = street->GetNumberInflowRate(b); // #/s
                          }
                        if (this->option_method == "ETR")
                          {
                            StreetNetworkTransport<T>::
                              ETRConcentration(transfer_velocity,
                                               temp,
                                               outgoing_flux,
                                               street_volume,
                                               number_concentration_array,
                                               number_concentration_array_tmp,
                                               background_number_concentration_array,
                                               number_emission_rate_array,
                                               number_inflow_rate_array,
                                               deposition_flow_rate_array_aer,
                                               street_deposition_flow_rate_array_aer,
                                               scavenging_flow_rate_array_aer,
                                               resuspension_factor,
                                               washoff_factor_array,
                                               street_surface_deposited_number_array,
                                               new_number_concentration_array,
                                               sub_delta_t_init,
                                               this->Nbin_aer,
                                               zero,
                                               new_street_surface_deposited_number_array);
                          }
                        else if(this->option_method == "Rosenbrock")
                          {
                            StreetNetworkTransport<T>::
                              RosenbrockConcentration(transfer_velocity,
                                                      temp,
                                                      outgoing_flux,
                                                      street_volume,
                                                      number_concentration_array,
                                                      number_concentration_array_tmp,
                                                      background_number_concentration_array,
                                                      number_emission_rate_array,
                                                      number_inflow_rate_array,
                                                      deposition_flow_rate_array_aer,
                                                      new_number_concentration_array,
                                                      sub_delta_t_init,
                                                      inflow_flux,
                                                      street_deposition_flow_rate_array_aer,
                                                      scavenging_flow_rate_array_aer,
                                                      resuspension_factor,
                                                      washoff_factor_array_aer,
                                                      street_surface_deposited_number_array,
                                                      new_street_surface_deposited_number_array,
                                                      this->Nbin_aer,
                                                      zero);
                          }
                        for (int b = 0; b < this->Nbin_aer; ++b)
                          {
                            street->SetStreetNumberConcentration(new_number_concentration_array(b), b);
                            street->SetStreetSurfaceDepositedNumber(new_street_surface_deposited_number_array(b), b);
                          }
                      }
                    if(number_computation_option == "based_on_mass")
                      {
                        for (int b = 0; b < this->Nbin_aer; ++b)
                          {
                            Data<T, 1> concentration_aer_bin(this->Ns_aer);
                            concentration_aer_bin.SetZero();
                            T TotalMass, Rho_aer;
                            T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um
                            TotalMass = 0.0;
                            for (int s = 0; s < this->Ns_aer - 1; s++)
                              {
                                concentration_aer_bin(s) = street->GetStreetConcentration_aer(s, b);
                                TotalMass += concentration_aer_bin(s);
                              }

                            // Fixed density or not
                            if (this->option_process["with_fixed_density"])
                              // conversion from kg/m3 to µg/µm3
                              Rho_aer = fixed_density_aer * 1.e-9;
                            else
                              Rho_aer = ComputeDensity(concentration_aer_bin, Rho_species, TotalMass, this->Ns_aer);
                            
                            new_number_concentration_array(b) = TotalMass/Rho_aer/pi*6./(MeanDiameter*MeanDiameter*MeanDiameter);

                            street->SetStreetNumberConcentration(new_number_concentration_array(b), b);
                          }
                      }
                  }

            
                //! Calculates the new sub_delta_t for the next iteration
                T sub_delta_t_max = next_date.GetSecondsFrom(next_date_tmp); 
                sub_delta_t_max = max(sub_delta_t_max, 0.0);

                if (this->option_process["with_number_concentration"] and
                    number_computation_option == "based_on_transport")
                  AdaptTimeStep(new_concentration_array,
                                concentration_array_tmp,
                                new_concentration_array_aer,
                                concentration_array_aer_tmp,
                                new_number_concentration_array,
                                number_concentration_array_tmp,
                                sub_delta_t_init,
                                this->sub_delta_t_min,
                                sub_delta_t_max,
                                sub_delta_t);
                else
                  AdaptTimeStep(new_concentration_array,
                                concentration_array_tmp,
                                new_concentration_array_aer,
                                concentration_array_aer_tmp,
                                sub_delta_t_init,
                                this->sub_delta_t_min,
                                sub_delta_t_max,
                                sub_delta_t);

              }

            
	    //! Chemical reactions
	    if (this->option_process["with_chemistry"])
	      {
		Array<T, 1> wet_diameter_aer(this->Nbin_aer);
		wet_diameter_aer = 0.0;
		for (int b = 0; b < this->Nbin_aer; ++b)
		  wet_diameter_aer(b) = street->GetStreetWetDiameter_aer(b);

		Array<T, 1> wet_diameter_aer_loc = wet_diameter_aer;
		
		Array<T, 1> photolysis_rate_array(this->Nr_photolysis);
		photolysis_rate_array = 0.0;

		for (int r = 0; r < this->Nr_photolysis; r++)
		  photolysis_rate_array(r) = street->GetPhotolysisRate(r);


		T attenuation_ = street->GetAttenuation();
		T specific_humidity_ = street->GetSpecificHumidity();
		T pressure_ = street->GetPressure();
		T temperature_ = street->GetTemperature();
		T longitude_ = street->GetLongitude();
		T latitude_ = street->GetLatitude();
		T rain_ = street->GetRain();
		T liquidwatercontent_ = street->GetLiquidWaterContent();
                
		Chemistry(current_date_tmp,
			  sub_delta_t_init,
			  attenuation_,
			  specific_humidity_,
			  pressure_,
			  temperature_,
			  longitude_,
			  latitude_,
			  rain_,
			  liquidwatercontent_,
			  photolysis_rate_array,
			  new_concentration_array,
			  new_concentration_array_aer,
			  new_number_concentration_array,
			  wet_diameter_aer,
			  wet_diameter_aer_loc);

		if (wet_diameter_option == "chemistry")
		  for (int b = 0; b < this->Nbin_aer; ++b)
		    street->SetStreetWetDiameter_aer(wet_diameter_aer_loc(b), b);
	      }

	    //! Set the new concentrations.
	    for (int s = 0; s < this->Ns; ++s)
	      street->SetStreetConcentration(new_concentration_array(s), s);

	    for (int b = 0; b < this->Nbin_aer; ++b)
	      {
		street->SetStreetNumberConcentration(new_number_concentration_array(b), b);
		for (int s = 0; s < this->Ns_aer; ++s)
		  {
		    street->SetStreetConcentration_aer(new_concentration_array_aer(s, b), s, b);
		  }
	      }

	    //! Update current_time_tmp
	    current_date_tmp.AddSeconds(sub_delta_t_init);
	    next_date_tmp.AddSeconds(sub_delta_t);
	    
	    //! Set the new sub_delta_t
	    sub_delta_t_init = sub_delta_t;
	  }

	for (int s = 0; s < this->Ns; ++s)
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
          {
            // The following arrays are updated using new computed data by each process.
            // They will be communicated by all processes later.
            concentration_mpi(i, s) = new_concentration_array(s);
	    massflux_roof_to_background_mpi(i, s) = temp * (new_concentration_array(s) - background_concentration_array(s)); // ug/s
            T conc_delta = new_concentration_array(s) - init_concentration_array(s);
            street_quantity_delta_mpi(i, s) = conc_delta * street_volume; // ug
          }
#else
	  {
            T massflux_roof_to_background;
	    massflux_roof_to_background = temp * (new_concentration_array(s) - background_concentration_array(s)); // ug/s
	    
	    street->SetMassfluxRoofToBackground(massflux_roof_to_background, s);
   
	    T conc_delta = new_concentration_array(s) - init_concentration_array(s);
	    T street_quantity_delta = conc_delta * street_volume; // ug
	    street->SetStreetQuantityDelta(street_quantity_delta, s);
	  }
#endif
          
          
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
            {
              // The following arrays are updated using new computed data by each process.
              // They will be communicated by all processes later.
              concentration_aer_mpi(i, s, b) = new_concentration_array_aer(s, b);
              massflux_roof_to_background_aer_mpi(i, s, b) = temp * (new_concentration_array_aer(s, b) - street->GetBackgroundConcentration_aer(s, b)); // ug/s
              T conc_delta = new_concentration_array_aer(s, b) - init_concentration_array_aer(s, b);
              street_quantity_delta_aer_mpi(i, s, b) = conc_delta * street_volume; // ug
            }
#else
	    {
	      T massflux_roof_to_background;
	      massflux_roof_to_background = temp * (new_concentration_array_aer(s, b) - street->GetBackgroundConcentration_aer(s, b)); // ug/s
	    
	      street->SetMassfluxRoofToBackground_aer(massflux_roof_to_background, s, b);
   
	      T conc_delta = new_concentration_array_aer(s, b) - init_concentration_array_aer(s, b);
	      T street_quantity_delta = conc_delta * street_volume; // ug
	      street->SetStreetQuantityDelta_aer(street_quantity_delta, s, b);
	    }
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
            for (int b = 0; b < this->Nbin_aer; ++b)
              concentration_number_mpi(i, b) = new_number_concentration_array(b);
#endif            
      }


#ifdef POLYPHEMUS_PARALLEL_WITH_MPI

    // Gather data from processes to Process 0
    // and cast all data from Process 0 to all other processes.
    BaseModuleParallel::GatherSlice_source_MPI(concentration_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(massflux_roof_to_background_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(street_quantity_delta_mpi);

    BaseModuleParallel::GatherSlice_source_MPI(concentration_aer_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(massflux_roof_to_background_aer_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(street_quantity_delta_aer_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(concentration_number_mpi);

    
    // Update data from MPI arrays to Street class objects.
    for (int i = 0; i < nstreet; i++)
      {
        Street<T>* street = this->StreetVector.at(i);
        for (int s = 0; s < this->Ns; s++)
          {
            street->SetStreetConcentration(concentration_mpi(i, s), s);
            street->SetMassfluxRoofToBackground(massflux_roof_to_background_mpi(i, s), s);
            street->SetStreetQuantityDelta(street_quantity_delta_mpi(i, s), s);
          }

        for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
            {
              street->SetStreetConcentration_aer(concentration_aer_mpi(i, s, b), s, b);
	      street->SetMassfluxRoofToBackground_aer(massflux_roof_to_background_aer_mpi(i, s, b), s, b);
              street->SetStreetQuantityDelta_aer(street_quantity_delta_aer_mpi(i, s, b), s, b);              
            }

        for (int b = 0; b < this->Nbin_aer; ++b)
          street->SetStreetNumberConcentration(concentration_number_mpi(i, b), b);

        
      }
#endif

  }

  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>
  ::AdaptTimeStep(const Array<T, 1> new_concentration,
		  const Array<T, 1>concentration_tmp,
		  const Array<T, 2> new_concentration_aer,
		  const Array<T, 2>concentration_tmp_aer,
		  const Array<T, 1> new_number_concentration,
		  const Array<T, 1>number_concentration_tmp,
		  const T sub_delta_t_init,
		  const T sub_delta_t_min,
		  const T sub_delta_t_max,
		  T& sub_delta_t)
  {
    T tmp, R;
    T EPSER = 0.01; //relative error precision
    T n2err = 0.0;
    
    //! local error estimation
    for (int s = 0; s < this->Ns; s++)
      if(new_concentration(s) > 0.0)
	{
	  tmp = (new_concentration(s) - concentration_tmp(s))/new_concentration(s);
	  n2err = n2err + tmp*tmp;
	}

    for (int b = 0; b < this->Nbin_aer; b++)
      {
	if(new_number_concentration(b) > 0.0)
	  {
	    tmp = (new_number_concentration(b) - number_concentration_tmp(b))/new_number_concentration(b);
	    n2err = n2err + tmp*tmp;
	  }
	for (int s = 0; s < this->Ns_aer - 1; s++) //do not consider H2O
	  if(new_concentration_aer(s,b) > 0.0)
	    {
	      tmp = (new_concentration_aer(s,b) - concentration_tmp_aer(s,b))/new_concentration_aer(s,b);
	      n2err = n2err + tmp*tmp;
	    }
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

  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>
  ::AdaptTimeStep(const Array<T, 1> new_concentration,
		  const Array<T, 1>concentration_tmp,
		  const Array<T, 2> new_concentration_aer,
		  const Array<T, 2>concentration_tmp_aer,
		  const T sub_delta_t_init,
		  const T sub_delta_t_min,
		  const T sub_delta_t_max,
		  T& sub_delta_t)
  {
    T tmp, R;
    T EPSER = 0.01; //relative error precision
    T n2err = 0.0;
    
    //local error estimation
    for (int s = 0; s < this->Ns; s++)
      if(new_concentration(s) > 0.0)
	{
	  tmp = (new_concentration(s) - concentration_tmp(s))/new_concentration(s);
	  n2err = n2err + tmp*tmp;
	}

    for (int s = 0; s < this->Ns_aer - 1; s++) //do not consider H2O
      for (int b = 0; b < this->Nbin_aer; b++)
	if(new_concentration_aer(s,b) > 0.0)
	  {
	    tmp = (new_concentration_aer(s,b) - concentration_tmp_aer(s,b))/new_concentration_aer(s,b);
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

  
  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::InitAllData()
  {
    StreetNetworkChemistry<T, ClassChemistry>::InitAllData();

    /*** Additional meteo data for aerosol chemistry on the streets ***/
    string filename;
    if (this->option_process["with_local_data"])
      {
	//  *** Meteo for the streets ***/
	this->InitData("meteo",
		       "LiquidWaterContent",
		       FileLiquidWaterContent_i,
		       FileLiquidWaterContent_f,
		       this->current_date,
		       LiquidWaterContent_f);
      }

    /*** Aerosol wet diameter ***/
    
    Array<T, 1> WetDiameter_aer(this->Nbin_aer);
    WetDiameter_aer = 0.0;
    InitWetDiameter_aer(WetDiameter_aer);

    /*** Particles density : unit g/cm3 ***/

    for (int s = 0; s < this->Ns_aer; s++)
      Rho_species_aer.push_back(ssh_mass_density_layers(s));
    
    /*** Emission data ***/
    int i, j, st;
    string species, species_bin;
    vector<int> isize_section;
    Date date;
    T Delta_t;
    date = this->input_files["emission_aer"].GetDateMin();
    Delta_t
      = this->input_files["emission_aer"].GetDelta_t();

    //! External mixing
    if (aerosol_emission_format == "Internal" && Nfraction_aer > 1)
      {
	Emission_aer_i.SetZero();
	Emission_aer_f.SetZero();
	FileEmission_aer_i.SetZero();
	FileEmission_aer_f.SetZero();
	NumberEmission_aer_i.SetZero();
	NumberEmission_aer_f.SetZero();
	FileNumberEmission_aer_i.SetZero();
	FileNumberEmission_aer_f.SetZero();
      }
    
    for (i = 0; i < Ns_emis_aer; i++)
      {
	species = species_list_emis_aer[i].first;
	isize_section= species_list_emis_aer[i].second;
	int CompositionID =
          FindExternalCompositionID(species);
	for (j = 0; j < int(isize_section.size()); j++)
	  {
	    species_bin = species + string("_") + to_str(isize_section[j]);
	    string filename
	      = this->input_files["emission_aer"](species_bin);

            //! External mixing using internal mixing format data.
	    if (aerosol_emission_format == "Internal" && Nfraction_aer > 1)
	      {
		this->InitData(filename, date, Delta_t,
			       FileEmission_aer_i,
			       FileEmission_aer_f,
			       this->current_date,
			       i,
			       j,
			       CompositionID,
			       Emission_aer_f);
	      }
	    else
	      {
		this->InitData(filename, date, Delta_t,
			       FileEmission_aer_i,
			       FileEmission_aer_f,
			       this->current_date,
			       i,
			       j,
			       Emission_aer_f,
			       this->Ncomposition_aer);		  
	      }   
	  }
      }

    //! Read or Compute aerosol emission data in number
    if (this->option_process["with_number_concentration"])
      for( j = 0; j < Nb_emis_aer; j++)
	{

          //! Read aerosol emission data in number if available.
          string filename = "";
          //! This option is set to false by default.
          //! If you have Number_X.bin files,
          //! Please add the following line in your data configuration
          //! With_emission_number_data: yes
          if (this->option_process["with_emission_number_data"])
            {
              species_bin = "Number_"  + to_str(emis_bin_list_aer[j]);
              filename = this->input_files["emission_aer"](species_bin);
            }
          
	  if (exists(filename))
	    {
	      if (aerosol_emission_format == "Internal" &&
                  Nfraction_aer > 1)
		{
		  TinyVector<int, 1> new_shape;
		  for (int i = 0; i < 1; i++)
		    new_shape(i) = NumberEmission_aer_f.GetArray().shape()(i + 1);

		  Data<T, 1> FileData_extract_i(new_shape);
		  Data<T, 1> FileData_extract_f(new_shape);
		  Data<T, 1> CurrentData_extract(new_shape);

		  this->InitData(filename, date, Delta_t,
				 FileData_extract_i,
				 FileData_extract_f,
				 this->current_date,
				 CurrentData_extract);
		  
		  for(st=0; st<this->total_nstreet; st++)
		    {
		      T TotalMass=0;
		      Data<T, 1> CompositionMass(this->Ncomposition_aer);
		      CompositionMass.SetZero();
		      for(int s=0; s<Ns_emis_aer; s++)
			{
			  species = species_list_emis_aer[s].first;
			  int CompositionID=FindExternalCompositionID(species);
			  int RealID=this->Ncomposition_aer*j+CompositionID;
			  TotalMass+=Emission_aer_f(s,RealID,st);
			  CompositionMass(CompositionID)+=Emission_aer_f(s,RealID,st);
			}
		      for(int id=0; id<this->Ncomposition_aer; id++)
			{
			  if(TotalMass>0)
			    {
			      int RealID=this->Ncomposition_aer*j+id;
			      NumberEmission_aer_f(RealID,st)=CurrentData_extract(st)*CompositionMass(id)/TotalMass;
			      FileNumberEmission_aer_i(RealID,st)=FileData_extract_i(st)*CompositionMass(id)/TotalMass;
			      FileNumberEmission_aer_f(RealID,st)=FileData_extract_f(st)*CompositionMass(id)/TotalMass;
			    }
			}
		    }
		}
	      else
		{
		  this->InitData(filename, date, Delta_t,
				 FileNumberEmission_aer_i,
				 FileNumberEmission_aer_f,
				 this->current_date, j,
				 NumberEmission_aer_f,
				 this->Ncomposition_aer);
		}
	    }
	  else //need to be calculate at last
	    {
	      ComputeNumberEmission_aer(int(emis_bin_list_aer[j]));
	    }
	}

    /*** Background concentration ***/
    
    if (this->option_process["with_local_data"])
      
      {
	date = this->input_files["bg_concentration_aer"].GetDateMin();
	Delta_t
	  = this->input_files["bg_concentration_aer"].GetDelta_t();
	if(aerosol_bg_format=="Internal"&&Nfraction_aer > 1)
	  {
	    Background_aer_i.SetZero();
	    Background_aer_f.SetZero();
	    FileBackground_aer_i.SetZero();
	    FileBackground_aer_f.SetZero();
	    NumberBackground_aer_i.SetZero();
	    NumberBackground_aer_f.SetZero();
	    FileNumberBackground_aer_i.SetZero();
	    FileNumberBackground_aer_f.SetZero();
	  }
	for (i = 0; i < Ns_bg_aer; i++)
	  {
	    species = species_list_bg_aer[i].first;
	    isize_section= species_list_bg_aer[i].second;
	    int CompositionID=FindExternalCompositionID(species);
	    for (j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["bg_concentration_aer"](species_bin);
		if(aerosol_bg_format=="Internal"&&Nfraction_aer > 1)
		  {
		    this->InitData(filename, date, Delta_t,
				   FileBackground_aer_i,
				   FileBackground_aer_f,
				   this->current_date,
				   i,
				   j,
				   CompositionID,
				   Background_aer_f);
		  }
		else
		  {
		    this->InitData(filename, date, Delta_t,
				   FileBackground_aer_i,
				   FileBackground_aer_f,
				   this->current_date,
				   i,
				   j,
				   Background_aer_f,
				   this->Ncomposition_aer);		  
		  }
	      }
	  }

        //! Get background number concentration.
	if (this->option_process["with_number_concentration"])
	  for( j = 0; j < Nb_bg_aer; j++)
	    {
	      species_bin = "Number_"  + to_str(bg_bin_list_aer[j]);
              string filename = this->input_files["bg_number_concentration"](species_bin);
	      if (exists(filename) and
                  this->option_process["with_bg_number_concentration_data"])
		{
		  if(aerosol_bg_format=="Internal"&&Nfraction_aer > 1)
		    {
		      TinyVector<int, 1> new_shape;
		      for (int i = 0; i < 1; i++)
			new_shape(i) = NumberBackground_aer_f.GetArray().shape()(i+1);
		      Data<T, 1> FileData_extract_i(new_shape);
		      Data<T, 1> FileData_extract_f(new_shape);
		      Data<T, 1> CurrentData_extract(new_shape);
		      this->InitData(filename, date, Delta_t,
				     FileData_extract_i,
				     FileData_extract_f,
				     this->current_date,
				     CurrentData_extract);
		      
		      for(st=0; st<this->total_nstreet; st++)
			{
			  T TotalMass=0;
			  Data<T, 1> CompositionMass(this->Ncomposition_aer);
			  CompositionMass.SetZero();
			  for(int s=0; s<Ns_bg_aer; s++)
			    { 
			      species = species_list_bg_aer[s].first;
			      int CompositionID=FindExternalCompositionID(species);
			      int RealID=this->Ncomposition_aer*j+CompositionID;
			      TotalMass+=Background_aer_f(s,RealID,st);
			      CompositionMass(CompositionID)+=Background_aer_f(s,RealID,st);
			    }
			  for(int id=0; id<this->Ncomposition_aer; id++)
			    {
			      if(TotalMass>0)
				{
				  int RealID=this->Ncomposition_aer*j+id;
				  NumberBackground_aer_f(RealID,st)=CurrentData_extract(st)*CompositionMass(id)/TotalMass;
				  FileNumberBackground_aer_i(RealID,st)=FileData_extract_i(st)*CompositionMass(id)/TotalMass;
				  FileNumberBackground_aer_f(RealID,st)=FileData_extract_f(st)*CompositionMass(id)/TotalMass;
				}
			    }
			}
		    }
		  else
		    { 
		      this->InitData(filename, date, Delta_t,
				     FileNumberBackground_aer_i,
				     FileNumberBackground_aer_f,
				     this->current_date, j,
				     NumberBackground_aer_f, this->Ncomposition_aer);
		    }
		}
	      else //need to be calculate at last
		ComputeNumberBackground_aer
		  (int(bg_bin_list_aer[j]));
	    }
      }// with local data
    
    /*** Road traffic ***/
    if (this->option_process["with_resuspension"])
      {
	this->InitData("traffic",
		       "RoadTraffic_2R",
		       FileRoadTraffic_2R_i,
		       FileRoadTraffic_2R_f,
		       this->current_date,
		       RoadTraffic_2R_f);
	this->InitData("traffic",
		       "RoadTraffic_HDV",
		       FileRoadTraffic_HDV_i,
		       FileRoadTraffic_HDV_f,
		       this->current_date,
		       RoadTraffic_HDV_f);
	this->InitData("traffic",
		       "RoadTraffic_PC",
		       FileRoadTraffic_PC_i,
		       FileRoadTraffic_PC_f,
		       this->current_date,
		       RoadTraffic_PC_f);
	this->InitData("traffic",
		       "RoadTraffic_LCV",
		       FileRoadTraffic_LCV_i,
		       FileRoadTraffic_LCV_f,
		       this->current_date,
		       RoadTraffic_LCV_f);
      }

    if (this->option_process["with_initial_condition_aer"])
      {
        string species, species_bin;
        vector<int> isize_section;

	StreetConcentration_aer.SetZero();
	if (ic_format == "Internal" && 
            this->option_process["with_external_composition"])
          {
            // necessary to translate into external composition
            StreetNumberConcentration_i.SetZero();
          }

        for (int i = 0; i < Ns_ic_aer; i++)
          {
            species = species_list_ic_aer[i].first;
            isize_section = species_list_ic_aer[i].second;

            for (int j = 0; j < int(isize_section.size()); j++)
              {
                species_bin = species + string("_") + to_str(isize_section[j]);
                string filename
                  = this->input_files["initial_condition_aer"](species_bin);
                int index = this->GetSpeciesIndex_aer(species);
                if (ic_format == "Internal" && 
                    this->option_process["with_external_composition"])
                  {
                    Data<T, 1> Concentration_tmp(&StreetConcentration_aer_i
                                                 (index, isize_section[j],0),
                                                 shape(this->total_nstreet));
                    if (is_num(filename))
                      Concentration_tmp.Fill(to_num<T>(filename));
                    else
                      FormatBinary<float>().Read(filename, Concentration_tmp);

                  }
		else
                  {
                    Data<T, 2> Concentration_tmp(&StreetConcentration_aer
                                                 (index,
                                                  isize_section[j] * this->Ncomposition_aer,
                                                  0),
                                                 shape(this->Ncomposition_aer, this->total_nstreet));
                    if (is_num(filename))
                      Concentration_tmp.Fill(to_num<T>(filename));
                    else
                      FormatBinary<float>().Read(filename, Concentration_tmp);

                  }
	      }
	  }
      }


    if (this->option_process["with_initial_condition_number_aer"])
      {
        StreetNumberConcentration.SetZero();
        for (int i = 0; i < Nb_ic_aer; i++)
          {
            species_bin = string("Number_") + to_str(ic_bin_list_aer[i]);
            string filename = this->input_files["initial_condition_aer"](species_bin);//problem of None
            if (exists(filename))
              {
                if(ic_format == "Internal" &&
                   this->option_process["with_external_composition"])
                  {
                    Data<T, 1> NumberConcentration_tmp(&StreetNumberConcentration_i
                                                       (ic_bin_list_aer[i], 0),
                                                       shape(this->total_nstreet));
                    if (is_num(filename))
                      NumberConcentration_tmp.Fill(to_num<T>(filename));
                    else 
                      FormatBinary<float>().Read(filename, NumberConcentration_tmp);
                  }
                else 
                  {
                    Data<T, 2> NumberConcentration_tmp(&StreetNumberConcentration
                                                       (ic_bin_list_aer[i] * this->Ncomposition_aer, 0),
                                                       shape(this->Ncomposition_aer, this->total_nstreet));
                    if (is_num(filename))
                      NumberConcentration_tmp.Fill(to_num<T>(filename));
                    else
                      FormatBinary<float>().Read(filename, NumberConcentration_tmp);
                  }
              }
          }
      }

    if (ic_format == "Internal" &&
        this->option_process["with_external_composition"])
      {
        // Set the initial condition to the first bin of each size.
        int Nfrac = this->Nbin_aer / this->Nsize_section_aer;
        for (int j = 0; j < this->Nsize_section_aer; j++)
          for(int ist = 0; ist < this->total_nstreet; ist++)
            for (int i = 0; i < Ns_ic_aer; i++)
              {
                species = species_list_ic_aer[i].first;
                int index_species = this->GetSpeciesIndex_aer(species);
                StreetConcentration_aer(index_species, j * Nfrac, ist) =
                  StreetConcentration_aer_i(index_species, j, ist);
                if (this->option_process["with_initial_condition_number_aer"])
                  {
                    StreetNumberConcentration(j * Nfrac, ist) =
                      StreetNumberConcentration_i(j, ist);
                  }                            
              }
      }
                 
    if (this->option_process["with_chemistry"])
      {
        Chemistry_.InitDistribution(this->StreetConcentration,
                                    StreetConcentration_aer,
                                    StreetNumberConcentration);
      }

    InitStreetConc();
    
  }


  //! Update the concentration arrays for each street.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::InitStreetConc()
  {

    StreetNetworkTransport<T>::InitStreetConc();

    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        
        for (int s = 0; s < this->Ns_aer; ++s)
          for (int b = 0; b < this->Nbin_aer; ++b)
              street->SetStreetConcentration_aer(StreetConcentration_aer(s, b, ist),
                                                 s, b);


        for (int b = 0; b < this->Nbin_aer; ++b)
          street->SetStreetNumberConcentration(StreetNumberConcentration(b, ist),
                                               b); 
          
        ++ist;
      }
  }  


  
//   ! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::InitStep()
  {
    StreetNetworkChemistry<T, ClassChemistry>::InitStep();
    int i, j, st, b;
    string species, species_bin;
    vector<int> isize_section;
    Date date;
    T Delta_t;

    //  *** Meteo for the streets ***/
    if (this->option_process["with_local_data"])
      this->UpdateData("meteo",
		       "LiquidWaterContent",
		       FileLiquidWaterContent_i,
		       FileLiquidWaterContent_f,
		       LiquidWaterContent_i,
		       LiquidWaterContent_f);

    /*** Aerosol wet diameter ***/
    Array<T, 1> WetDiameter_aer(this->Nbin_aer);
    WetDiameter_aer = 0.0;
    if (wet_diameter_option == "gerber")
      InitWetDiameter_aer(WetDiameter_aer);

    /*** Aerosol emissions ***/
    date = this->input_files["emission_aer"].GetDateMin();
    Delta_t
      = this->input_files["emission_aer"].GetDelta_t();	
 
    for (i = 0; i < Ns_emis_aer; i++)
      {
	species = species_list_emis_aer[i].first;
	isize_section = species_list_emis_aer[i].second;
	int CompositionID = FindExternalCompositionID(species);
	for (j = 0; j < int(isize_section.size()); j++)
	  {
	    species_bin = species + string("_") + to_str(isize_section[j]);
	    string filename
	      = this->input_files["emission_aer"](species_bin);
	    if(aerosol_emission_format == "Internal" && Nfraction_aer > 1)
	      {		  
		this->UpdateData(filename, date, Delta_t,
				 FileEmission_aer_i,
				 FileEmission_aer_f,
				 i, j,CompositionID,
				 Emission_aer_i,
				 Emission_aer_f);
	      }
	    else
	      {
		this->UpdateData(filename, date, Delta_t,
				 FileEmission_aer_i,
				 FileEmission_aer_f,
				 i,
				 j,
				 Emission_aer_i,
				 Emission_aer_f,
				 this->Ncomposition_aer);
	      }
	  }
      }  


    //! Read or Compute aerosol emission data in number    
    if (this->option_process["with_number_concentration"])
      for( j = 0; j < Nb_emis_aer ; j++)
	{

          //! Read aerosol emission data in number if available.
          string filename = "";
          //! This option is set to false by default.
          //! If you have Number_X.bin files,
          //! Please add the following line in your data configuration
          //! With_emission_number_data: yes
          if (this->option_process["with_emission_number_data"])
            {
              species_bin = "Number_"  + to_str(emis_bin_list_aer[j]);
              filename = this->input_files["emission_aer"](species_bin);
            }
          
	  if (exists(filename))
	    {
	      if (aerosol_emission_format == "Internal" &&
                  Nfraction_aer > 1)
		{
		  TinyVector<int, 1> new_shape;
		  for (i = 0; i < 1; i++)
		    new_shape(i) = NumberEmission_aer_f.GetArray().shape()(i + 1);

		  Data<T, 1> FileData_extract_i(new_shape);
		  Data<T, 1> FileData_extract_f(new_shape);
		  Data<T, 1> CurrentData_extract_i(new_shape);
		  Data<T, 1> CurrentData_extract_f(new_shape);
		  this->UpdateData(filename, date, Delta_t,
				   FileData_extract_i,
				   FileData_extract_f,
				   CurrentData_extract_i,
				   CurrentData_extract_f);
		  for(int st = 0; st < this->total_nstreet; st++)
		    {
		      T TotalMass = 0;
		      Data<T, 1> CompositionMass(this->Ncomposition_aer);
		      CompositionMass.SetZero();
		      for(int s = 0; s < Ns_emis_aer; s++)
			{
			  species = species_list_emis_aer[s].first;
			  int CompositionID = FindExternalCompositionID(species);
			  int RealID = this->Ncomposition_aer * j + CompositionID;
			  TotalMass += Emission_aer_f(s, RealID, st);
			  CompositionMass(CompositionID) += Emission_aer_f(s, RealID, st);
			}
		      for(int id = 0; id < this->Ncomposition_aer; id++)
			{
			  if(TotalMass>0)
			    {
			      int RealID = this->Ncomposition_aer * j + id;
			      NumberEmission_aer_i(RealID, st) =
                                CurrentData_extract_i(st) * CompositionMass(id) / TotalMass;
			      NumberEmission_aer_f(RealID, st) =
                                CurrentData_extract_f(st) * CompositionMass(id) / TotalMass;
			      FileNumberEmission_aer_i(RealID, st) =
                                FileData_extract_i(st) * CompositionMass(id) / TotalMass;
			      FileNumberEmission_aer_f(RealID, st) =
                                FileData_extract_f(st) * CompositionMass(id) / TotalMass;
			    }
			}
		    }
		}
	      else
		{
		  this->UpdateData(filename, date, Delta_t,
				   FileNumberEmission_aer_i,
				   FileNumberEmission_aer_f,
				   j,
				   NumberEmission_aer_i,
				   NumberEmission_aer_f,
				   this->Ncomposition_aer);
		}
	    }
	  else
	    {
	      for(int id = 0 ; id < this->Ncomposition_aer; id++ )
		for (st = 0; st < this->total_nstreet; st++)
		  NumberEmission_aer_i(j * this->Ncomposition_aer + id, st) =
		    NumberEmission_aer_f(j * this->Ncomposition_aer + id, st);
		    
	      ComputeNumberEmission_aer(int(emis_bin_list_aer[j]));
	    }
	}
   
    if (this->option_process["with_local_data"])
      { 
	/*** Aerosol background ***/
	date = this->input_files["bg_concentration_aer"].GetDateMin();
	Delta_t
	  = this->input_files["bg_concentration_aer"].GetDelta_t();	
	for (i = 0; i < Ns_bg_aer; i++)
	  {
	    species = species_list_bg_aer[i].first;
	    isize_section = species_list_bg_aer[i].second;
	    int CompositionID = FindExternalCompositionID(species);
	    for (j = 0; j < int(isize_section.size()); j++)
	      {
		species_bin = species + string("_") + to_str(isize_section[j]);
		string filename
		  = this->input_files["bg_concentration_aer"](species_bin);
		if (aerosol_bg_format == "Internal" && Nfraction_aer > 1)
		  {		  
		    this->UpdateData(filename, date, Delta_t,
				     FileBackground_aer_i,
				     FileBackground_aer_f,
				     i, j,CompositionID,
				     Background_aer_i,
				     Background_aer_f);
		  }
		else
		  {

		    this->UpdateData(filename, date, Delta_t,
				     FileBackground_aer_i,
				     FileBackground_aer_f,
				     i,
				     j,
				     Background_aer_i,
				     Background_aer_f,
				     this->Ncomposition_aer);
		  }
	      }
	  }

	if (this->option_process["with_number_concentration"])
	  for (j = 0; j < Nb_bg_aer ; j++)
	    {
	      species_bin = "Number_"  + to_str(bg_bin_list_aer[j]);
	      string filename = this->input_files["bg_number_concentration"](species_bin);
	      if (exists(filename))
		{
		  if (aerosol_bg_format == "Internal" && Nfraction_aer > 1)
		    {
		      TinyVector<int, 1> new_shape;
		      for (i = 0; i < 1; i++)
			new_shape(i) = NumberBackground_aer_f.GetArray().shape()(i + 1);

		      Data<T, 1> FileData_extract_i(new_shape);
		      Data<T, 1> FileData_extract_f(new_shape);
		      Data<T, 1> CurrentData_extract_i(new_shape);
		      Data<T, 1> CurrentData_extract_f(new_shape);
		      
		      this->UpdateData(filename, date, Delta_t,
				       FileData_extract_i,
				       FileData_extract_f,
				       CurrentData_extract_i,
				       CurrentData_extract_f);
		      
		      for (int st=0; st < this->total_nstreet; st++)
			{
			  T TotalMass=0;
			  Data<T, 1> CompositionMass(this->Ncomposition_aer);
			  CompositionMass.SetZero();
			  for (int s=0; s<Ns_bg_aer; s++)
			    {
			      species = species_list_bg_aer[s].first;
			      int CompositionID = FindExternalCompositionID(species);
			      int RealID = this->Ncomposition_aer * j + CompositionID;
			      TotalMass += Background_aer_f(s, RealID, st);
			      CompositionMass(CompositionID) +=
                                Background_aer_f(s, RealID, st);
			    }
			  for (int id = 0; id < this->Ncomposition_aer; id++)
			    {
			      if (TotalMass > 0)
				{
				  int RealID = this->Ncomposition_aer * j + id;
				  NumberBackground_aer_i(RealID, st) =
                                    CurrentData_extract_i(st) * CompositionMass(id) / TotalMass;
				  NumberBackground_aer_f(RealID, st) =
                                    CurrentData_extract_f(st) * CompositionMass(id) / TotalMass;
				  FileNumberBackground_aer_i(RealID, st) =
                                    FileData_extract_i(st) * CompositionMass(id) / TotalMass;
				  FileNumberBackground_aer_f(RealID, st) =
                                    FileData_extract_f(st) * CompositionMass(id) / TotalMass;
				}
			    }
			}
		    }
		  else
		    {

		      this->UpdateData(filename, date, Delta_t,
				       FileNumberBackground_aer_i,
				       FileNumberBackground_aer_f, j,
				       NumberBackground_aer_i,
				       NumberBackground_aer_f, this->Ncomposition_aer);
		    }
		}
	      else
		{
		  for(int id =0 ; id < this->Ncomposition_aer; id++ )
		    for (st = 0; st < this->total_nstreet; st++)
		      NumberBackground_aer_i(j * this->Ncomposition_aer + id, st) =
			NumberBackground_aer_f(j * this->Ncomposition_aer + id, st);
			
		  ComputeNumberBackground_aer(int(bg_bin_list_aer[j]));
		}
	    } 
      }

    /*** Road Traffic ***/
    if (this->option_process["with_resuspension"])
      {
	this->UpdateData("traffic",
			 "RoadTraffic_2R",
			 FileRoadTraffic_2R_i,
			 FileRoadTraffic_2R_f,
			 RoadTraffic_2R_i,
			 RoadTraffic_2R_f);
	this->UpdateData("traffic",
			 "RoadTraffic_HDV",
			 FileRoadTraffic_HDV_i,
			 FileRoadTraffic_HDV_f,
			 RoadTraffic_HDV_i,
			 RoadTraffic_HDV_f);
	this->UpdateData("traffic",
			 "RoadTraffic_PC",
			 FileRoadTraffic_PC_i,
			 FileRoadTraffic_PC_f,
			 RoadTraffic_PC_i,
			 RoadTraffic_PC_f);
	this->UpdateData("traffic",
			 "RoadTraffic_LCV",
			 FileRoadTraffic_LCV_i,
			 FileRoadTraffic_LCV_f,
			 RoadTraffic_LCV_i,
			 RoadTraffic_LCV_f);
	
	SetStreetRoadTraffic(RoadTraffic_2R_f,
			     RoadTraffic_HDV_f,
			     RoadTraffic_PC_f,
			     RoadTraffic_LCV_f);
      }

    Array<T, 2> emission_rate_aer(this->Ns_aer, this->Nbin_aer);
    Array<T, 1> number_emission_rate(this->Nbin_aer);    

    st = 0;   

    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;
	
        //! Set the emission data
        emission_rate_aer = 0.0; 
        for (i = 0; i < this->Ns_aer; ++i)
          for (j = 0; j < Ns_emis_aer; ++j)
            {
              if (this->species_list_aer[i] == species_list_emis_aer[j].first)
                {
                  for (b = 0; b < Nb_emis_aer; ++b)
                    {
                      emission_rate_aer(i, b) = Emission_aer_f(j, b, st);
                    }
                }
            }

        for (b = 0; b < Nb_emis_aer; ++b)
	  number_emission_rate(b) = NumberEmission_aer_f(b, st);

	street->SetEmission_aer(emission_rate_aer);
	street->SetNumberEmission(number_emission_rate);
	
        if (this->option_process["with_local_data"])
	  {	  
            for (i = 0; i < this->Ns_aer; ++i)
	      for (j = 0; j < Ns_bg_aer; ++j)
                if (this->species_list_aer[i] == species_list_bg_aer[j].first)
                  for (b = 0; b < Nb_bg_aer; ++b)
                    {
                      street->SetBackgroundConcentration_aer
                        (Background_aer_f(j, b, st), i, b);		      
                      if (!this->option_process["with_initial_condition_aer"])
                        {
                          //Concentration at streets equal to Cbg in first time step
                          if(this->current_date == this->Date_min)
                            street->SetStreetConcentration_aer
                              (Background_aer_f(j, b, st), i, b);
                        }
                    }
	    if (this->option_process["with_number_concentration"])
              for (b = 0; b < Nb_bg_aer; ++b)
                street->SetBackgroundNumberConcentration(NumberBackground_aer_f(b, st), b);

            //! Set the meteo data
	    street->SetLiquidWaterContent(LiquidWaterContent_f(st));
	  }	 	
	++st;	  
      }

    //! Initilialize the inflow rate.
    InitInflowRate();
  }


  //! Chemistry.
  /*!
    \param 
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::Chemistry()
  {
    Array<T, 1> source(this->Ns);
    Array<T, 1> photolysis_rate(this->Nr_photolysis);
    Array<T, 1> concentration(this->Ns);
    
    source = 0.0;

    int nfoglay = 0;
    T lwc_avg = 0;
    T heightfog = 0;
    int ifog = 0;
    Array<T, 1> incloudwetdepositionflux(this->Ns);
    Array<T, 2> incloudwetdepositionflux_aer(this->Ns_aer, this->Nbin_aer);
    Array<T, 1> incloudwetdepositionfluxnumber(this->Nbin_aer);
    Array<T, 2> concentration_aer(this->Ns_aer, this->Nbin_aer);
    incloudwetdepositionflux = 0.0;
    incloudwetdepositionflux_aer = 0.0;
    incloudwetdepositionfluxnumber = 0.0;
    Array<T, 1> number_concentration(this->Nbin_aer);
    int ninterface = 2;
    Array<T, 1> VerticalInterface(ninterface);
    VerticalInterface = 0.0;

    int st = 0;

    // MPI implementation.
    // 'For' loop on the street segments is parallelized.
    int first_index_along_source, last_index_along_source;
    
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI

    int nstreet = this->StreetVector.size();

    // The following arrays need to be used only for MPI.
    Array<T, 2> concentration_mpi;
    Array<T, 3> concentration_aer_mpi;
    Array<T, 2> concentration_number_mpi;
    Array<T, 2> wet_diameter_mpi;
    concentration_mpi.resize(nstreet, this->Ns);
    concentration_aer_mpi.resize(nstreet, this->Ns_aer, this->Nbin_aer);
    concentration_number_mpi.resize(nstreet, this->Nbin_aer);
    wet_diameter_mpi.resize(nstreet, this->Nbin_aer);
    
    // Scatter buffers for Street-type arrays to all processes.
    BaseModuleParallel::ScatterSlice_source_MPI(concentration_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(concentration_aer_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(concentration_number_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(wet_diameter_mpi);
    
    // Get the first and last indices of streets for each process.
    // The number of slices is the number of streets / the number of processes.
    // For example, if 1000 streets are computed by 10 processes,
    // Process 1: Street 1 (first_index) to 100 (last_index),
    // Process 2: Street 101 (first_index) to 200 (last_index), ...
    // See BaseModuleParallel.cxx
    BaseModuleParallel::GetEdgePartition_source
      (first_index_along_source, last_index_along_source);
#else
    first_index_along_source = 0;
    last_index_along_source = this->StreetVector.size();
#endif

    // Parallelized loop.
    for (int i = first_index_along_source; i < last_index_along_source; i++)
      {

        T pH = 4.5; 
        Street<T>* street = this->StreetVector.at(i);
	for (int s = 0; s < this->Ns; ++s)
          concentration(s) = street->GetStreetConcentration(s);
        
	for (int s = 0; s < this->Ns_aer; ++s)
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    concentration_aer(s,b) = street->GetStreetConcentration_aer(s, b);
       
	for (int b = 0; b < this->Nbin_aer; ++b)
	  {
	    number_concentration(b) = street->GetStreetNumberConcentration(b); 
	  }
        
        T attenuation_ = street->GetAttenuation();
        T specific_humidity_ = street->GetSpecificHumidity();
        T pressure_ = street->GetPressure();
        T temperature_ = street->GetTemperature();
        T longitude_ = street->GetLongitude();
        T latitude_ = street->GetLatitude();
	T rain_ = street->GetRain();
	T liquidwatercontent_ = street->GetLiquidWaterContent();

        for (int r = 0; r < this->Nr_photolysis; r++)
	  photolysis_rate(r) = street->GetPhotolysisRate(r);
	// ! Get the wet diameter
	Array<T, 1> wet_diameter_aer(this->Nbin_aer);
	wet_diameter_aer = 0.0;
	for (int b = 0; b < this->Nbin_aer; ++b)
	  wet_diameter_aer(b) = street->GetStreetWetDiameter_aer(b);

        // Call Forward in Aerosol_SSH.cxx
        Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
        		   attenuation_,
			   specific_humidity_,
        		   temperature_,
			   pressure_,
			   source,
        		   photolysis_rate,
        		   T(this->next_date.
        		     GetSecondsFrom(this->current_date)),
        		   attenuation_,
			   specific_humidity_,
        		   temperature_,
			   pressure_,
			   source,
        		   photolysis_rate,
			   longitude_, latitude_,
			   concentration,
			   liquidwatercontent_, 
			   wet_diameter_aer,
			   concentration_aer,
			   pH,
			   number_concentration);

        // Call Forward_aer in Aerosol_SSH.cxx
	Chemistry_.Forward_aer(T(this->current_date.GetNumberOfSeconds()), 
			       specific_humidity_,
			       temperature_,
			       pressure_,
			       T(this->next_date.GetSecondsFrom(this->current_date)),
			       concentration, 
			       liquidwatercontent_,
			       rain_,
			       VerticalInterface, 
			       concentration_aer,
			       incloudwetdepositionflux, 
			       incloudwetdepositionflux_aer,
			       pH,
			       lwc_avg,
			       heightfog,
			       ifog, 
			       number_concentration,
			       incloudwetdepositionfluxnumber,
                               wet_diameter_aer);


        /*!  The following arrays are updated using new computed data by each process.
          They will be communicated by all processes later.
        */
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
        {
          for (int s = 0; s < this->Ns; ++s)
            concentration_mpi(i, s) = concentration(s);

          for (int s = 0; s < this->Ns_aer; ++s)
            for (int b = 0; b < this->Nbin_aer; ++b)
              concentration_aer_mpi(i, s, b) = concentration_aer(s, b);
        
          for (int b = 0; b < this->Nbin_aer; ++b)
            concentration_number_mpi(i, b) = number_concentration(b);

          if (wet_diameter_option == "chemistry")
            for (int b = 0; b < this->Nbin_aer; ++b)
              wet_diameter_mpi(i, b) = wet_diameter_aer(b); 

        }
#else
        {
          for (int s = 0; s < this->Ns; ++s)
            street->SetStreetConcentration(concentration(s), s);
          
          for (int s = 0; s < this->Ns_aer; ++s)
            for (int b = 0; b < this->Nbin_aer; ++b)
              street->SetStreetConcentration_aer(concentration_aer(s, b), s, b); 

          for (int b = 0; b < this->Nbin_aer; ++b)
            street->SetStreetNumberConcentration(number_concentration(b), b); 

	if (wet_diameter_option == "chemistry")
	  for (int b = 0; b < this->Nbin_aer; ++b)
	    street->SetStreetWetDiameter_aer(wet_diameter_aer(b), b);
        }
#endif
	st += 1;
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI

    // Gather data from processes to Process 0
    // and cast all data from Process 0 to all other processes.
    BaseModuleParallel::GatherSlice_source_MPI(concentration_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(concentration_aer_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(concentration_number_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(wet_diameter_mpi);
    
    // Update data from MPI arrays to Street class objects.
    for (int i = 0; i < nstreet; i++)
      {
        Street<T>* street = this->StreetVector.at(i);
        for (int s = 0; s < this->Ns; s++)
          street->SetStreetConcentration(concentration_mpi(i, s), s);

        for (int s = 0; s < this->Ns_aer; ++s)
          for (int b = 0; b < this->Nbin_aer; ++b)
            street->SetStreetConcentration_aer(concentration_aer_mpi(i, s, b),
                                               s, b); 

        for (int b = 0; b < this->Nbin_aer; ++b)
          street->SetStreetNumberConcentration(concentration_number_mpi(i, b),
                                               b);

        for (int b = 0; b < this->Nbin_aer; ++b)
          street->SetStreetWetDiameter_aer(wet_diameter_mpi(i, b),
                                           b);

      }
#endif
  }
  
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::Chemistry(Date current_date_tmp,
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
							  Array<T, 1>& concentration,
							  Array<T, 2>& concentration_aer,
							  Array<T, 1>& number_concentration,
							  Array<T, 1> wet_diameter_aer,
							  Array<T, 1>& wet_diameter_aer_loc)
  {  

    Array<T, 1> source(this->Ns);
    source = 0.0;

    int nfoglay = 0;
    T pH = 4.5; 
    T lwc_avg = 0;
    T heightfog = 0;
    int ifog = 0;
    Array<T, 1> incloudwetdepositionflux(this->Ns);
    Array<T, 2> incloudwetdepositionflux_aer(this->Ns_aer, this->Nbin_aer);
    Array<T, 1> incloudwetdepositionfluxnumber(this->Nbin_aer);
    incloudwetdepositionflux = 0.0;
    incloudwetdepositionflux_aer = 0.0;
    incloudwetdepositionfluxnumber = 0.0;
    int ninterface = 2;
    Array<T, 1> VerticalInterface(ninterface);
    VerticalInterface = 0.0;
    
    Chemistry_.Forward(T(current_date_tmp.GetNumberOfSeconds()),
		       attenuation_,
		       specific_humidity_,
		       temperature_,
		       pressure_,
		       source,
		       photolysis_rate,
		       sub_delta_t,
		       attenuation_,
		       specific_humidity_,
		       temperature_,
		       pressure_,
		       source,
		       photolysis_rate,
		       longitude_, latitude_,
		       concentration,
		       liquidwatercontent_, 
		       wet_diameter_aer,
		       concentration_aer,
		       pH,
		       number_concentration);

    
    Chemistry_.Forward_aer(T(current_date_tmp.GetNumberOfSeconds()), 
			   specific_humidity_,
			   temperature_,
			   pressure_,
			   sub_delta_t,
			   concentration, 
			   liquidwatercontent_,
			   rain_,
			   VerticalInterface, 
			   concentration_aer,
			   incloudwetdepositionflux, 
			   incloudwetdepositionflux_aer,
			   pH,
			   lwc_avg,
			   heightfog,
			   ifog, 
			   number_concentration,
			   incloudwetdepositionfluxnumber,
    			   wet_diameter_aer_loc);
    
  }

  //in the first time step the wet diameter is set equal to dry diameter. The correct wet diameter is calculated in the chemical module
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::InitWetDiameter_aer
  (Array<T, 1>& WetDiameter_aer_)
  {
    int b, id;
    Array<T, 1> MeanDiameter(this->Nsize_section_aer);
    for (b = 0; b < this->Nsize_section_aer; b++)
      MeanDiameter(b) = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b));
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
	 iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;
	T relative_humidity_ = street->GetRelativeHumidity();
	T temperature_ = street->GetTemperature();
	for (b = 0; b < this->Nsize_section_aer; b++)
	  for (id = 0; id < this->Ncomposition_aer; id ++)
	    {    
	      int bin_id=b*this->Ncomposition_aer+id;
	      if(relative_humidity_>0.0)
		_gerber_wet_diameter(&relative_humidity_,
				     &temperature_, &MeanDiameter(b),
				     &WetDiameter_aer_(bin_id));
	      else
		WetDiameter_aer_(bin_id)=MeanDiameter(b);
	    
	      // Back to meters.
	      WetDiameter_aer_(bin_id) *= 1.e-6;
	      street->SetStreetWetDiameter_aer(WetDiameter_aer_(bin_id), bin_id);
	    }
      }
  }

  //! Sets the road traffic for the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetRoadTraffic(Data<T, 1> RoadTraffic_2R_f,
								     Data<T, 1> RoadTraffic_HDV_f,
								     Data<T, 1> RoadTraffic_PC_f,
								     Data<T, 1> RoadTraffic_LCV_f)
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;	  
	street->SetStreetRoadTraffic(RoadTraffic_2R_f(ist),
				     RoadTraffic_HDV_f(ist),
				     RoadTraffic_PC_f(ist),
				     RoadTraffic_LCV_f(ist));
        ++ist;
      }
  }
  
  //! Sets the aerosol concentrations for the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetConcentration_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int s = 0; s < this->Ns_aer; ++s)
	    for (int b = 0; b < this->Nbin_aer; ++b)
              StreetConcentration_aer(s, b, ist) = street->GetStreetConcentration_aer(s, b);

        ++ist;

      }
  }

   //! Sets the number concentrations for the whole street-network.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::SetStreetNumberConcentration_aer()
  {
    int ist = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        for (int b = 0; b < this->Nbin_aer; ++b)
          StreetNumberConcentration(b, ist) = street->GetStreetNumberConcentration(b);
        ++ist;
      }
  }  
  
  //! Compute number concentration from mass concentration, diameter and
  //! density
  /*!
    \ param b bin number
    \ return density (ug.um-3)
  */
  template<class T, class ClassChemistry>
  T StreetNetworkAerosol<T, ClassChemistry>
  ::ComputeDensity(Data <T, 1> Conc_aer_tmp,
                   vector<T> Rho_species, T TotalMass, int Ns)
  {
    int s;
    T rho, subrho;
    
    rho = 0.0;
    subrho = 0.0;
	
    for (s = 0; s < Ns; s++)
      subrho += Conc_aer_tmp(s) / Rho_species[s];

    if (TotalMass < 1.e-10 or subrho == 0.)
      rho = 1.e-6;
    else
      rho = 1.e-6 * TotalMass/subrho;

    if (rho == 0.0)
      throw string("Error: zero density in ComputeDenisty.\n") +
        "Total concentration is " + to_str(TotalMass) +
        " and sum of the ratio of conc/density is " + to_str(subrho) + "." +
        "Plese check if the density is set to zero for a species.";
    
    return rho;	
  }     
  
  
  //! Compute number volume emission from mass concentration, diameter and
  //! density
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>
  ::ComputeNumberEmission_aer(int b)
  {
    //here b is the index of size section
    //id_b is the index of bin
    int id_b = this->NumberEmissionIndex_aer(b*this->Ncomposition_aer);
    int index_b = Bin_to_size_index_aer(id_b);	
    int st;
    T TotalMass; // ug.m-3
    T Rho_aer; 
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um
    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species; // g.cm-3 -> 10-6 ug.um-3
	
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc=this->Ncomposition_aer;
    
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(this-> Ns_emis_aer);
    Conc_aer_tmp.SetZero();
	
    for (int s = 0; s < Ns_emis_aer; s++)
      {
        rho_stream_aer.PeekValue(species_list_emis_aer[s].first, Rho_tmp);
        Rho_species.push_back(Rho_tmp);
      }

    for (st = 0; st < this->total_nstreet; st++)
      {
        for (int ic = 0; ic < Nc; ic++) //ZS
          {
            TotalMass = 0.0;
            Rho_aer = 0.0;
            for (int s = 0; s < Ns_emis_aer; s++)
              {
                it_begin = species_list_emis_aer[s].second.begin();
                it_end = species_list_emis_aer[s].second.end();
                pos = find(it_begin,it_end,b);
                dist = distance(it_begin, pos);

                if (dist < int(species_list_emis_aer[s].second.size()))
                  {
                    TotalMass = TotalMass +
                      this->Emission_aer_f(s,dist*Nc+ic,st);
                    Conc_aer_tmp(s) = this->Emission_aer_f(s,dist*Nc+ic,st);
                  }

              }

            // Fixed density or not
            if (this->option_process["with_fixed_density"])
              // conversion from kg/m3 to µg/µm3
              Rho_aer = fixed_density_aer * 1.e-9;
            else
              Rho_aer = ComputeDensity(Conc_aer_tmp,
                                       Rho_species, TotalMass, Ns_emis_aer);
            
            this->NumberEmission_aer_f(index_b*Nc+ic,st)=
              TotalMass/Rho_aer/pi*6.
              /(MeanDiameter*MeanDiameter*MeanDiameter); 
            double tmp_n=this->NumberEmission_aer_f(index_b*Nc+ic,st);
            if(TotalMass*tmp_n==0&&TotalMass!=tmp_n)
              cout<<st<<" , "<<" m:"<<TotalMass<<" n:"<<tmp_n
                  <<" Rho:"<<Rho_aer<<" MeanDiameter:"<<MeanDiameter<<endl;
	    
          }
      }
  }
  
    //! Compute number volume emission from mass concentration, diameter and
  //! density
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>
  ::ComputeNumberBackground_aer(int b)
  {
    //here b is the index of size section
    //id_b is the index of bin
    int id_b = this->NumberBackgroundIndex_aer(b*this->Ncomposition_aer);
    int index_b = Bin_to_size_index_aer(id_b);	
    int st;
    T TotalMass; // ug.m-3
    T Rho_aer; 
    T MeanDiameter = 1.e6 * sqrt(BinBound_aer(b + 1) * BinBound_aer(b)); // um
    string saver_file;
    T Rho_tmp;
    vector<T> Rho_species; // g.cm-3 -> 10-6 ug.um-3
	
    this->config.SetSection("[data]");
    this->config.PeekValue("Data_description", saver_file);
    ConfigStream saver_stream_Nb(saver_file);
    ConfigStream rho_stream_aer(this->file_species);
    rho_stream_aer.SetSection("[mass_density_aer]");
    vector<int>::iterator pos;
    vector<int>::iterator it_begin;
    vector<int>::iterator it_end;
    int dist;
    int Nc=this->Ncomposition_aer;
    
    Data<T, 1> Conc_aer_tmp;
    Conc_aer_tmp.Resize(Ns_bg_aer);
    Conc_aer_tmp.SetZero();
	
    for (int s = 0; s < Ns_bg_aer; s++)
      {
	rho_stream_aer.PeekValue(species_list_bg_aer[s].first, Rho_tmp);
	Rho_species.push_back(Rho_tmp);
      }
	
    for (st = 0; st < this->total_nstreet; st++)
      {
	for (int ic = 0; ic < Nc; ic++) //ZS
	{
	  TotalMass = 0.0;
	  Rho_aer = 0.0;
	  for (int s = 0; s < Ns_bg_aer; s++)
	    {
	      it_begin = species_list_bg_aer[s].second.begin();
	      it_end = species_list_bg_aer[s].second.end();
	      pos = find(it_begin,it_end,b);
	      dist = distance(it_begin, pos);

	      if (dist < int(species_list_bg_aer[s].second.size()))
		{
		  TotalMass = TotalMass + Background_aer_f(s,dist*Nc+ic,st);
		  Conc_aer_tmp(s) = Background_aer_f(s,dist*Nc+ic,st);
		}

	    }

          // Fixed density or not
          if (this->option_process["with_fixed_density"])
            // conversion from kg/m3 to µg/µm3
            Rho_aer = fixed_density_aer * 1.e-9;
          else
            Rho_aer = ComputeDensity(Conc_aer_tmp, Rho_species, TotalMass, Ns_bg_aer);
          
	  this->NumberBackground_aer_f(index_b*Nc+ic,st)=
	    TotalMass/Rho_aer/pi*6.
	    /(MeanDiameter*MeanDiameter*MeanDiameter);
	  double tmp_n=this->NumberBackground_aer_f(index_b*Nc+ic,st);	    
	}
      }
  }  
  
  //! Write the results in the output files in the text format.
  template<class T, class ClassChemistry>
  void StreetNetworkAerosol<T, ClassChemistry>::OutputSaver()
  {
    StreetNetworkTransport<T>::OutputSaver();
  }
  
  /////////////////
  // External composition methods //SZ
  /////////////////

  //! Return external composition id based on species name
  //! (where this specie related group is pure)
  /*!
    \ param b bin number
    \ return number concentration (m-3)
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>::FindExternalCompositionID(string species)
  {
    int index_species = this->GetSpeciesIndex_aer(species);
    int index_group = this->aerosol_species_group_relation(index_species);
    
    if(index_group==(this->Ngroup_aer-1))
      return 0;
    else
      for(int i = 1; i < this->Ncomposition_aer; i++)
      {
	if(this->composition_bounds(i,index_group,1)==1.0)
	  return i;
      }
    cout<<"Error!:can not find related composition ID,";
    cout<<"please check species-group relations"<<endl;
    return 0;
  }  
  
  //! Checks whether the model deals with number concentration
  /*!
    \param s species global index.
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template<class T, class ClassChemistry>
  bool StreetNetworkAerosol<T, ClassChemistry>::HasNumberConcentration_aer()
  {
    if(this->option_process["with_number_concentration"])
      return true;
    else
      return false;
  }  
  
  //! Returns the concentrations for the whole street-network.
  template<class T, class ClassChemistry>
  inline Data<T, 3>& StreetNetworkAerosol<T, ClassChemistry>::GetStreetConcentration_aer()
  {
    return StreetConcentration_aer;
  }  
  
  //! Returns the concentrations for the whole street-network.
  template<class T, class ClassChemistry>
  inline Data<T, 2>& StreetNetworkAerosol<T, ClassChemistry>::GetStreetNumberConcentration()
  {
    return StreetNumberConcentration;
  }  
  
  //! Returns the index of size section based on the index of bin .
  /*!
    \param b bin number.
    \return index of size section.
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>
  ::Bin_to_size_index_aer(int b) const
  {
    if(b<this->Nbin_aer && b >=0)
      return (b/this->Ncomposition_aer);
    else
      throw string("Error: Bin index") + string(":") + to_str(b)
	+ "is not exit!";
  }  
  
  //! Returns the index in number volume emissions of a given aerosol.
  /*!
    \param b bin number.
    \return The aerosol index in volume emissions.
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>
  ::NumberEmissionIndex_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    if(HasNumberEmission_aer(x))
    {
      int binout=find(emis_bin_list_aer.begin(),
		  emis_bin_list_aer.end(), b)
	- emis_bin_list_aer.begin();
      return Bin_index_translate_aer(binout,x);
    }
    else
      throw string("Species \"Number") + string("_") + to_str(b)
	+ "\" not found in Volume Emissions.";
  }  
  
  //! Checks whether an aerosol has number volume emissions.
  /*!
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template<class T, class ClassChemistry>
  bool StreetNetworkAerosol<T, ClassChemistry>
  ::HasNumberEmission_aer(int x) const
  {//b is size section index, x is bin index
    int b=Bin_to_size_index_aer(x);
    return find(this->emis_bin_list_aer.begin(),
		this->emis_bin_list_aer.end(), b)
      != this->emis_bin_list_aer.end();
  }  
  
  //! translate bin b from current size section to new size section s.
  /*!
    \param b bin number. s the new size section number
    \return index of compsition section.
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>
  ::Bin_index_translate_aer(int s, int b) const
  {
    int composition_id=Bin_to_composition_index(b);
    if(b<this->Nbin_aer && b >=0)
      return (s*this->Ncomposition_aer+composition_id);
    else
      throw string("Error: Bin index") + string(":") + to_str(b)
	+ "is not exit!";
  }  
  
  //! Returns the index of compsition section based on the index of bin .
  /*!
    \param b bin number.
    \return index of compsition section.
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>
  ::Bin_to_composition_index(int b) const
  {
    if(b<this->Nbin_aer && b >=0)
      return (b-Bin_to_size_index_aer(b)*this->Ncomposition_aer);
    else
      throw string("Error: Bin index") + string(":") + to_str(b)
	+ "is not exit!";
  }  
  
  //! Returns the index in number volume emissions of a given aerosol.
  /*!
    \param b bin number.
    \return The aerosol index in volume emissions.
  */
  template<class T, class ClassChemistry>
  int StreetNetworkAerosol<T, ClassChemistry>
  ::NumberBackgroundIndex_aer(int x) const
  {
    int b=Bin_to_size_index_aer(x);
    if(HasNumberBackground_aer(x))
    {
      int binout=find(bg_bin_list_aer.begin(),
		  bg_bin_list_aer.end(), b)
	- bg_bin_list_aer.begin();
      return Bin_index_translate_aer(binout,x);
    }
    else
      throw string("Species \"Number") + string("_") + to_str(b)
	+ "\" not found in Volume Emissions.";
  }   
  
  //! Checks whether an aerosol has number volume emissions.
  /*!
    \param b bin number.
    \return True if the aerosol has volume emissions, false otherwise.
  */
  template<class T, class ClassChemistry>
  bool StreetNetworkAerosol<T, ClassChemistry>
  ::HasNumberBackground_aer(int x) const
  {//b is size section index, x is bin index
    int b=Bin_to_size_index_aer(x);
    return find(this->bg_bin_list_aer.begin(),
		this->bg_bin_list_aer.end(), b)
      != this->bg_bin_list_aer.end();
  }    

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKAEROSOL_CXX
#endif
