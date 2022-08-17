// Copyright (C) 2020, CEREA (ENPC - EDF R&D)
// Author(s): Youngseob Kim
//
// This file is a component of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
//
// Polyphemus is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// Polyphemus is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the Polyphemus web site:
//      http://cerea.enpc.fr/polyphemus/


#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SSH_CXX

#include "Aerosol_SSH.hxx"

#define PI 3.141592653589793115997963468544185161590576171875
// Perfect gas constant ([atm.L/mol/K]).
#define Pr 8.2e-2
// Pi / 6.0
#define cst_pi6 0.52359877559829887307

namespace Polyphemus
{

  //! Default constructor.
  template<class T>
  Aerosol_SSH<T>::Aerosol_SSH()
  {
  }

  //! Destructor.
  template<class T>
  Aerosol_SSH<T>::~Aerosol_SSH()
  {

    //! SSH-aerosol shared library is closed
    //! ONLY if it has been open.
    if (_aerosol_so != NULL)
      {
        api.call(_aerosol_so, "api_sshaerosol_finalize_");

        char *error = NULL;

        // Close the shared library
        dlclose(_aerosol_so);
        error = dlerror();
        if (error != NULL)
          cout << "Error on dlclose :\t" << error << endl;
      }
  }
  

  //! Read ssh-aerosol shared library.
  template<class T>
  template<class ClassModel>
  void Aerosol_SSH<T>::InitSharedLib(ClassModel& Model)
  {
    // This is used to handle errors
    char *error = NULL;
    dlerror();    /* Clear any existing error */

    char lib_path[] = "./libssh-aerosol.so";
    
    // Try to load the shared library
    _aerosol_so = dlopen(lib_path, RTLD_LAZY);
    error = dlerror();
    if (error != NULL)
      cout << "Error on dlopen :\t" << error << endl;

    // Declare SSH-aerosol as not running standalone
    // This removes most of the output to stdout
    api.send_bool(_aerosol_so,
                  "api_sshaerosol_set_standalone_",
                  false);

    // Declare SSH-aerosol logs informations.
    api.send_bool(_aerosol_so,
                  "api_sshaerosol_set_logger_",
                  true);


    // Initialize SSH-aerosol with namelist.ssh
    char namelist_ssh[] = "namelist.ssh";
    api.exchange_char_array(_aerosol_so,
                            "api_sshaerosol_initialize_",
                            namelist_ssh);

    // Initialize the photolysis
    api.call(_aerosol_so, "api_sshaerosol_initphoto_");
    
    // Get the number of gas species
    Ns = api.recv_int(_aerosol_so,
                      "api_sshaerosol_get_ngas_");
    
    // Get the number of aerosol species
    Ns_aer_nolayer = api.recv_int(_aerosol_so,
                                  "api_sshaerosol_get_naero_");
    
    // Get the number of aerosol layers
    N_layer = api.recv_int(_aerosol_so,
                           "api_sshaerosol_get_nlayer_");

    // Get i_hydrophilic
    i_hydrophilic = api.recv_int(_aerosol_so,
                                 "api_sshaerosol_get_i_hydrophilic_");
    
    // Get the number of aerosol species in different layers
    Ns_aer = api.recv_int(_aerosol_so,
                          "api_sshaerosol_get_n_aerosol_layers_");

    // Get the number of aerosol bins
    Nbin_aer = api.recv_int(_aerosol_so,
                            "api_sshaerosol_get_nsize_");
    
    // Get the number of aerosol size-sections
    Nsize_section_aer = api.recv_int(_aerosol_so,
                                     "api_sshaerosol_get_nsizebin_");
    ssh_diam_input.resize(Nsize_section_aer + 1);
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_diam_input_",
                              ssh_diam_input);
    ssh_diam_input *= 1.e-6; // conversion from microm to m.

    // Get the number of photolysis
    Nr_photolysis = api.recv_int(_aerosol_so,
                                 "api_sshaerosol_get_nphotolysis_");

    // Whether the gas-phase chemistry taken into accout.
    with_gas_chemistry = api.recv_bool(_aerosol_so,
                                       "api_sshaerosol_get_tag_chem_");

    // Whether externally mixed.
    with_external_composition = api.recv_bool(_aerosol_so,
                                              "api_sshaerosol_get_tag_external_");

    // in g/mol
    ssh_mass_density_layers.resize(Ns_aer);
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_mass_density_layers_",
                              ssh_mass_density_layers);

    // Get the number of inorganic aerosol species
    Ns_inorganic_aer = api.recv_int(_aerosol_so,
                                    "api_sshaerosol_get_nesp_isorropia_");

    // Get the number of organic aerosol species
    Ns_organic_aer = api.recv_int(_aerosol_so,
                                  "api_sshaerosol_get_nesp_aec_");

    // Get the number of inert aerosol species
    Ns_inert_aer = api.recv_int(_aerosol_so,
                                "api_sshaerosol_get_n_inert_");

    // Get index of aerosol groups.
    index_groups.resize(Ns_aer_nolayer);
    api.exchange_int_array(_aerosol_so,
                           "api_sshaerosol_get_index_groups_",
                           index_groups);
 
    index_groups_ext.resize(this->Ns_aer);
    api.exchange_int_array(_aerosol_so,
                           "api_sshaerosol_get_index_groups_ext_",
                           index_groups_ext);    
    

    Ncomposition_aer = api.recv_int(_aerosol_so,
                                    "api_sshaerosol_get_ncomposition_");


    Nfraction_aer = api.recv_int(_aerosol_so,
                                 "api_sshaerosol_get_nfrac_");    

    Ngroup_aer = api.recv_int(_aerosol_so,
                              "api_sshaerosol_get_n_groups_");

    discretization_composition.resize(2, Ngroup_aer, Ncomposition_aer); 
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_discretization_composition_",
                              discretization_composition);

    aerosol_type.resize(Ns_aer_nolayer);
    api.exchange_int_array(_aerosol_so,
                           "api_sshaerosol_get_aerosol_type_",
                           aerosol_type);
    
    // Get names of aerosol species.
    aerosol_spec_name.resize(Ns_aer_nolayer);
    for (int n = 0; n < Ns_aer_nolayer; n++)
      {
        char aero_name[81];
        api.exchange_char_array(_aerosol_so,
                                "api_sshaerosol_get_aero_name_",
                                n + 1, aero_name);
        aerosol_spec_name(n) = trim(aero_name);
      }
    
    
  }

  
  //! Check configuration.
  template<class T>
  template<class ClassModel>
  void Aerosol_SSH<T>::CheckConfiguration(ClassModel& Model)
  {

    if (Model.GetNs() != Ns)
      throw string("Incompatibility: model manages ") + to_str(Model.GetNs())
	+ string(" species while chemistry has ") + to_str(Ns) + " species.";
    // if (Model.GetNr_photolysis() != Nr_photolysis)
    //   throw string("Incompatibility: model manages ")
    //     + to_str(Model.GetNr_photolysis())
    //     + string(" photolysis reactions while chemistry has ")
    //     + to_str(Nr_photolysis) + " photolysis reactions.";
    if (Model.GetNs_aer() != Ns_aer)
      throw string("Incompatibility: model manages ") + to_str(Model.GetNs_aer())
	+ string(" aerosol species while chemistry has ")
        + to_str(Ns_aer) + " species.";
    if (Model.GetNsize_section_aer() != Nsize_section_aer)
      throw string("Incompatibility: model manages ") + to_str(Model.GetNsize_section_aer())
	+ string(" size sections while chemistry has ") + to_str(Nsize_section_aer) + " size sections.";
     
  
  }

  
  //! Initialization of the scheme.
  /*! In case the model is incompatible with the chemical mechanism, an
    exception is thrown.
    \param Model model with the following interface:
    <ul>
    <li> GetNs()
    <li> GetNbin_aer()
    <li> GetNgroup_aer()
    <li> GetNcomposition_aer()
    <li> GetNsize_section_aer()
    <li> GetSpeciesList()
    <li> GetSpeciesFile()
    <li> GetNr_photolysis()
    <li> GetPhotolysisReactionList()
    <li> GetNs_source()
    <li> GetNz_source()
    <li> SourceGlobalIndex(int)
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SSH<T>::Init(ClassModel& Model)
  {

    ConfigStream config_species(Model.GetSpeciesFile());

    string species;

    // Determine gas species particle interacting index ().
    config_species.SetSection("[gas_species_cloud_interact]");

    while (!config_species.IsEmpty())
      {
	species = config_species.GetElement();
	int gas_index = Model.GetSpeciesIndex(species) + 1;
	species_index_cloud_interact.push_back(gas_index);
      }
    Ns_cloud_interact = int(species_index_cloud_interact.size());
    
    // Source splitting for volumic sources.
    Ns_source = Model.GetNs_source();
    Nz_source = Model.GetNz_source();
    source_index.resize(Ns_source);
    for (int i = 0; i < Ns_source; i++)
      source_index(i) = Model.SourceGlobalIndex(i);

    // Module options.
    ConfigStream config(Model.GetConfigurationFile());
    config.SetSection("[options]");

    if (config.Check("With_in_cloud_scavenging"))
      config.PeekValue
        ("With_in_cloud_scavenging",
         option_process_aer["with_in_cloud_scavenging"]);
    else
      // Default option is false
      option_process_aer["with_in_cloud_scavenging"] = false;  
    if (config.Check("Collect_wet_flux_aerosol"))
        config.PeekValue("Collect_wet_flux_aerosol",
                         option_process_aer["collect_wet_flux_aer"]);
    else
      // Default option is false      
      option_process_aer["collect_wet_flux_aer"] = false;
    if (config.Check("Collect_wet_flux"))
      config.PeekValue("Collect_wet_flux",
                       option_process_aer["collect_wet_flux"]);
    else
      // Default option is false      
      option_process_aer["collect_wet_flux"] = false;
    config.PeekValue("aqueous_module",aqueous_module);
    aqueous_module = lower_case(aqueous_module);
    config.PeekValue("With_number_concentration",
		     option_process_aer["with_number_concentration"]);
    // if (!option_process_aer["with_number_concentration"])
    //   throw string("Error: please set With_number_concentration to yes. ") +
    //     "This option is indispensable for SSH-aerosol module."; // YK
    config.PeekValue("With_fixed_density",
		     option_process_aer["with_fixed_density"]);
  
    // Density in kg / m^3.
    config.PeekValue("Fixed_aerosol_density", FixedDensity_aer);
    
    // Lwc cloud threshold.
    if (config.Check("Lwc_cloud_threshold"))
      config.PeekValue("Lwc_cloud_threshold", lwc_cloud_threshold);
    else
      lwc_cloud_threshold = 0.05;

    // Redistribution method for aqueous chemistry module
    config.PeekValue("Redistribution_method", redistribution_method);
    redistribution_method = lower_case(redistribution_method);

    // 
    if (config.Check("Conserving_mass_tolerance"))
      config.PeekValue("Conserving_mass_tolerance", "positive",
                       conserving_mass_tolerance);
    else
      conserving_mass_tolerance = 0.0;
    
    // Reads bin bounds.
    vector<string> bin_list;
    config.SetSection("[domain]");
    config.Find("Bin_bounds");
    bin_list = split(config.GetLine());

    // Reads bin bounds in micrometers and converts it to meters.
    BinBound_aer.resize(Model.GetNsize_section_aer() + 1);
    for (int i = 0; i < Model.GetNsize_section_aer() + 1; i++)
      BinBound_aer(i) = 1.e-6 * convert<T>(bin_list[i]);


    for (int i = 0; i < Nsize_section_aer + 1; i++)
      {
        if (ssh_diam_input(i) != BinBound_aer(i))
          throw string("Initial aerosol size is not same between the model ") + to_str(BinBound_aer(i))
            + string(" and SSH-aerosol ") + to_str(ssh_diam_input(i));
      }
 
    

    if (redistribution_method == "moving-diameter")
      iredist = 10;
    else if (redistribution_method == "siream")
      iredist = 11;
    else if (redistribution_method == "siream-euler-coupled")
      iredist = 12;
    else if (redistribution_method == "siream-moving-diameter")
      iredist = 13;
    else
      throw string("bad string \"") + redistribution_method +
	string("\" for redistribution method option,\n") +
	string("possibilities are moving-diameter, siream, siream-euler-coupled,");
   
    if (option_process_aer["with_in_cloud_scavenging"])
      {
	if (option_process_aer["collect_wet_flux_aer"])
	  {
	    if (!Model.HasField("InCloudWetDepositionFlux_aer"))
	      throw string("if field \" InCloudWetDepositionFlux_aer \" ") +
		string("does not exist in Model, in cloud wet deposition ") +
		string("fluxes cannot be collected.");
	    if (Model.D4("InCloudWetDepositionFlux_aer").GetLength(2)
		!= Model.GetNy()
		|| Model.D4("InCloudWetDepositionFlux_aer").GetLength(3)
		!= Model.GetNx())
	      throw string("\" InCloudWetDepositionFlux_aer \" has") +
		string(" wrong dimension.");
	    if (option_process_aer["with_number_concentration"])
	      {
		if(!Model.HasField("InCloudWetDepositionFluxNumber_aer"))
		  throw string("if field \" InCloudWetDepositionFluxNumber_aer \" ") +
		    string("does not exist in Model, in cloud wet deposition ") +
		    string(" number fluxes cannot be collected.");
		if (Model.D3("InCloudWetDepositionFluxNumber_aer").GetLength(1)
		    != Model.GetNy()
		    || Model.D3("InCloudWetDepositionFluxNumber_aer").GetLength(2)
		    != Model.GetNx())
		  throw string("\" InCloudWetDepositionFluxNumber_aer \" has") +
		    string(" wrong dimension.");
	      }
			
	  }
	if (option_process_aer["collect_wet_flux"])
	  {
	    if (!Model.HasField("InCloudWetDepositionFlux"))
	      throw string("if field \" InCloudWetDepositionFlux \" ") +
		string("does not exist in Model, in cloud wet deposition ") +
		string("fluxes cannot be collected.");
	    if (Model.D3("InCloudWetDepositionFlux").GetLength(1)
		!= Model.GetNy()
		|| Model.D3("InCloudWetDepositionFlux").GetLength(2)
		!= Model.GetNx())
	      throw string("\" InCloudWetDepositionFlux \" has") +
		string(" wrong dimension.");
	  }

	// Check if in-cloud scavenging takes place
	// and if an aqueous-phase module is not activated.
	// In-cloud scavenging must be accompanied with
        // an aqueous-phase module (simple or VSRM module).
     	if (aqueous_module=="no")
	  throw string("Warning: activate an aqueous module ") +
	    string("(simple or VSRM) in your configuration file, ") +
	    string("polair3d.cfg. ") +
	    string("In-cloud scavenging is calculated ") +
            string("in the aqueous module. ") +
	    string("See aerosol.f and Aerosol_SCRAM.cxx");
      }
 
    BaseModuleParallel::Init(Model);
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Partition along axis x.
    BuildPartition_x();
#endif

    /*** Display options ***/
    config.SetSection("[display]");
    // Should the configuration be displayed on screen?
    config.PeekValue("Show_configuration",
		     option_process_aer["show_configuration"]);
    
    if (option_process_aer["show_configuration"] && GetRank() == 0)
      DisplayConfiguration(Model);

    if (GetRank() == 0)
      {
        CheckConfiguration(Model);
      }
  }
  

  //! Display the configuration.
  template<class T>
  template<class ClassModel>
  void Aerosol_SSH<T>::DisplayConfiguration(ClassModel& Model)
  {

    ConfigStream config_species(Model.GetSpeciesFile());

    string species;

    // Determine gas species particle interacting index ().
    config_species.SetSection("[gas_species_cloud_interact]");
    cout << "Interacting cloud species : ";
    while (!config_species.IsEmpty())
      {
        species = config_species.GetElement();
        int gas_index = Model.GetSpeciesIndex(species) + 1;
        cout << species << "(" << gas_index << ") ";
      }
    cout << endl;

    cout << "Module: Fixed_aerosol_density = " << FixedDensity_aer << endl;
    cout << "Module: Lwc_cloud_threshold = " << lwc_cloud_threshold << endl;
    cout << "Module: redistribution method for aqueous chemistry: " <<
      redistribution_method << endl;
    
    cout << "\t==== SSH-aerosol module configuration ==== " << endl;
    cout << "SSH: number of gas-phase species = "<< Ns << endl;    
    cout << "SSH: number of photolysis = " << Nr_photolysis << endl;    
    cout << "SSH: number of aerosol species = " << Ns_aer << endl;
    cout << "SSH: number of aerosol layers = " << N_layer << endl;
    cout << "SSH: number of aerosol species without considering layers = " << Ns_aer_nolayer << endl;
    cout << "SSH: number of aerosol size section = " << Nsize_section_aer << endl;
    cout << "SSH: number of size-composition aerosol bin = " << Nbin_aer << endl;
    if (with_gas_chemistry)
      cout << "SSH: Run with gas-phase chemistry." << endl;
    else
      cout << "SSH: Run without gas-phase chemistry." << endl;
    
  }

    
  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("Attenuation_i")
    <li> D3("SpecificHumidity_i")
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetSource_i()
    <li> D4("PhotolysisRate_i")
    <li> D3("Attenuation_f")
    <li> D3("SpecificHumidity_f")
    <li> D3("Temperature_f")
    <li> D3("Pressure_f")
    <li> GetSource_f()
    <li> D4("PhotolysisRate_f")
    <li> GetGridXArray1D()
    <li> GetGridYArray1D()
    <li> GetConcentration()
    <li> D3("LiquidWaterContent_i")
    <li> D4("WetDiameter_aer")
    <li> GetConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SSH<T>::Forward(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();


    Forward(T(date_i.GetNumberOfSeconds()), Model.D3("Attenuation_i"),
            Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
            Model.D3("Pressure_i"), Model.GetSource_i(),
            Model.D4("PhotolysisRate_i"), T(date_f.GetSecondsFrom(date_i)),
            Model.D3("Attenuation_f"), Model.D3("SpecificHumidity_f"),
            Model.D3("Temperature_f"), Model.D3("Pressure_f"),
            Model.GetSource_f(), Model.D4("PhotolysisRate_f"),
            Model.GetGridXArray1D(), Model.GetGridYArray1D(),
            Model.GetConcentration(), Model.D3("LiquidWaterContent_i"),
            Model.D4("WetDiameter_aer"), Model.GetConcentration_aer(),
            Model.D3("pH"), Model.GetNumberConcentration_aer());
	
  }


  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("SpecificHumidity_i")
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetConcentration()
    <li> D3("LiquidWaterContent_i")
    <li> D2("Rain_i")
    <li> GetLayerInterface()
    <li> GetConcentration_aer()
    <li> D3("pH")
    <li> GetNumberConcentration_aer()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Aerosol_SSH<T>::Forward_aer(ClassModel& Model)
  {

    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();
	
    Data<T, 3> in_cloud_wet_flux(Ns, Model.GetNy(), Model.GetNx());
    Data<T, 4> in_cloud_wet_flux_aer(Ns_aer, Nbin_aer, Model.GetNy(),
				     Model.GetNx());
    Data<T, 3> in_cloud_wet_flux_number_aer(Nbin_aer, Model.GetNy(), Model.GetNx());

    Forward_aer(T(date_i.GetNumberOfSeconds()),
        	Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
        	Model.D3("Pressure_i"),	T(date_f.GetSecondsFrom(date_i)),
        	Model.GetConcentration(), Model.D3("LiquidWaterContent_i"),
        	Model.D2("Rain_i"), Model.GetLayerInterface(),
        	Model.GetConcentration_aer(), in_cloud_wet_flux,
        	in_cloud_wet_flux_aer,
        	Model.D3("pH"),
        	Model.GetNumberConcentration_aer(),
        	in_cloud_wet_flux_number_aer);
	
    if (option_process_aer["with_in_cloud_scavenging"])
      {
	if (option_process_aer["collect_wet_flux"])
	  for (int s = 0; s < Ns; s++)
	    if(Model.HasScavenging(s))
	      for (int j = 0; j < Model.GetNy(); j++)
		for (int i = 0; i < Model.GetNx(); i++)
		  Model.D3("InCloudWetDepositionFlux")
		    (Model.ScavengingIndex(s), j, i) =
		    in_cloud_wet_flux(s, j, i);
	if (option_process_aer["collect_wet_flux_aer"])
	  for (int b = 0; b < Nbin_aer; b++)
	    if(Model.HasScavenging_aer(b))
	      for (int j = 0; j < Model.GetNy(); j++)
		for (int i = 0; i < Model.GetNx(); i++)
		  {
		    for (int s = 0; s < Ns_aer; s++)
		      Model.D4("InCloudWetDepositionFlux_aer")
			(s, Model.ScavengingIndex_aer(b), j, i) =
			in_cloud_wet_flux_aer(s, b, j, i);
		    if(option_process_aer["with_number_concentration"])
		      Model.D3("InCloudWetDepositionFluxNumber_aer")
			(Model.ScavengingIndex_aer(b), j, i) =
			in_cloud_wet_flux_number_aer(b, j, i);
		  }
		  
      }
		
  }


  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param Attenuation_i cloud attenuation coefficients at the beginning of
    the time step.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param Source_i volume sources at the beginning of the time step.
    \param PhotolysisRate_i photolysis rates at the beginning of the time
    step.
    \param next_time time at the end of the time step.
    \param Attenuation_f cloud attenuation coefficients at the end of the
    time step.
    \param SpecificHumidity_f specific humidity at the end of the time step.
    \param Temperature_f temperature at the end of the time step.
    \param Pressure_f pressure at the end of the time step.
    \param Source_f volume sources at the end of the time step.
    \param PhotolysisRate_f photolysis rates at the end of the time step.
    \param Longitude longitudes.
    \param Latitude latitudes.
    \param Concentration gas concentrations.
    \param LiquidWaterContent_i air liquid water content at the beginning of
    the time step.
    \param WetDiameter_aer Aerosol wet diameters at the beginning of the time
    step.
    \param Concentration_aer aerosol concentrations.
  */
  template<class T>
  void Aerosol_SSH<T>::Forward(T current_time,
                               Data<T, 3>& Attenuation_i,
                               Data<T, 3>& SpecificHumidity_i,
                               Data<T, 3>& Temperature_i,
                               Data<T, 3>& Pressure_i,
                               Data<T, 4>& Source_i,
                               Data<T, 4>& PhotolysisRate_i,
                               T delta_t,
                               Data<T, 3>& Attenuation_f,
                               Data<T, 3>& SpecificHumidity_f,
                               Data<T, 3>& Temperature_f,
                               Data<T, 3>& Pressure_f,
                               Data<T, 4>& Source_f,
                               Data<T, 4>& PhotolysisRate_f,
                               Array<T, 1>& Longitude,
                               Array<T, 1>& Latitude,
                               Data<T, 4>& Concentration,
                               Data<T, 3>& LiquidWaterContent_i,
                               Data<T, 4>& WetDiameter_aer,
                               Data<T, 5>& Concentration_aer,
                               Data<T, 3>& pH,
                               Data<T, 4>& NumberConcentration_aer)
  
  {
   
    int iheter = option_process_aer["with_heterogeneous_reactions"] ? 1 : 0;

    int first_index_along_x, last_index_along_x;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // This array will store the concentration data with the species
    // axis permuted with the X-axis.
    Array<T,4> OrdConc;
    Array<T,5> OrdConc_aer;
    Array<T,4> OrdNumConc_aer;
    ScatterSlice_x_MPI(Concentration, OrdConc);
    ScatterSlice_x_MPI(Concentration_aer, OrdConc_aer);
    if (option_process_aer["with_number_concentration"])
      ScatterSlice_x_MPI(NumberConcentration_aer, OrdNumConc_aer);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
#else
    first_index_along_x = 0;
    last_index_along_x = Concentration.GetLength(3);
#endif
    
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)	\
  firstprivate(first_index_along_x, last_index_along_x)
#endif
    for (int i = first_index_along_x; i < last_index_along_x; i++)
      for (int k = 0; k < Concentration.GetLength(1); k++)
        for (int j = 0; j < Concentration.GetLength(2); j++)
	  {
            // Make output directories to save ssh-aerosol results.
            // It can be used only for one grid cell for the test purpose.

            Array<T, 1> Photolysis1D(Nr_photolysis);
            Array<T, 1> Photolysis1D_f(Nr_photolysis);
            Array<T, 1> Concentration1D(Ns);
            Array<T, 1> Source1D(Ns_source);
            Array<T, 1> Source1D_f(Ns_source);
            Array<T, 1> WetDiameter_aer1D(Nbin_aer);
            Array<T, 2> Concentration_aer2D(Ns_aer, Nbin_aer);
	    Array<T, 1> NumberConcentration_aer1D(Nbin_aer);
            int s, b;
           
	    for (s = 0; s < Nr_photolysis; s++)
	      {
		Photolysis1D(s) = PhotolysisRate_i()(s, k, j, i);
		Photolysis1D_f(s) = PhotolysisRate_f()(s, k, j, i);
	      }

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		Concentration1D(s) = OrdConc(i, k, j, s);
#else
		Concentration1D(s) = Concentration(s, k, j, i);
#endif
	      }

            for (b = 0; b < Nbin_aer ; b++)
              {
                WetDiameter_aer1D(b) = WetDiameter_aer(b, k, j, i);
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    Concentration_aer2D(s, b) = OrdConc_aer(i, b, k, j, s);
#else
                    Concentration_aer2D(s, b) = Concentration_aer(s, b, k, j, i);
#endif

                  }
              }

	    if (option_process_aer["with_number_concentration"])
	      for (b = 0; b < Nbin_aer; b++)
		{
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		  // Rank of indices must match the one defined in
		  // ScatterSlice_x_MPI.

		  if (isnan(OrdNumConc_aer(i, k, j, b)))
		    {
		      cout <<"isnan "<< b <<" "<< i <<" "<<  j<<" " <<  k;
						
		    }

					
		  NumberConcentration_aer1D(b) = OrdNumConc_aer(i, k, j, b);
#else

		  NumberConcentration_aer1D(b) = NumberConcentration_aer(b, k, j, i);
#endif
		}
	    else
	      NumberConcentration_aer1D = 0.0;

	    if (k < Nz_source)
	      {
		for (s = 0; s < Ns_source ; s++)
		  {
		    Source1D(s) = Source_i()(s, k, j, i);
		    Source1D_f(s) = Source_f()(s, k, j, i);
		  }
	      }
	    else
	      {
		Source1D = 0;
		Source1D_f = 0;
	      }
            
	    Forward(current_time, Attenuation_i(k, j, i),
		    SpecificHumidity_i(k, j, i), Temperature_i(k, j, i),
		    Pressure_i(k, j, i), Source1D, Photolysis1D,
		    delta_t, Attenuation_f(k, j, i),
		    SpecificHumidity_f(k, j, i), Temperature_f(k, j, i),
		    Pressure_f(k, j, i), Source1D_f, Photolysis1D_f,
		    Longitude(i), Latitude(j), Concentration1D, 
                    LiquidWaterContent_i(k, j, i), WetDiameter_aer1D, 
                    Concentration_aer2D, pH(k, j, i),
		    NumberConcentration_aer1D);

		    
	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		OrdConc(i, k, j, s) = Concentration1D(s);
#else
		Concentration(s, k, j, i) = Concentration1D(s);
#endif
	      }

	    for (s = 0; s < Ns_aer; s++)
              {
                for (b = 0; b < Nbin_aer; b++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    OrdConc_aer(i, b, k, j, s) = Concentration_aer2D(s, b);
#else
                    Concentration_aer(s, b, k, j, i) = Concentration_aer2D(s, b);
#endif
                  }
              }
	    if (option_process_aer["with_number_concentration"])
	      for (b = 0; b < Nbin_aer ; b++)
		{
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		  // Rank of indices must match the one defined in
		  // ScatterSlice_x_MPI.
		  OrdNumConc_aer(i, k, j, b) = NumberConcentration_aer1D(b);
#else
		  NumberConcentration_aer(b, k, j, i) = NumberConcentration_aer1D(b);
#endif
		} 


	  } // loop j
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(OrdConc, Concentration);
    GatherSlice_x_MPI(OrdConc_aer, Concentration_aer);
    if(option_process_aer["with_number_concentration"])
      GatherSlice_x_MPI(OrdNumConc_aer, NumberConcentration_aer);
#endif

  } // Forward


  
  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param next_time time at the end of the time step.
    \param Concentration gas concentrations.
    \param LiquidWaterContent_i air liquid water content at the beginning of
    the time step.
    \param Rain_i rain rate at the beginning of the time step.
    \param VerticalInterface height of vertical interfaces in meter.
    \param Concentration_aer aerosol concentrations.
    \param aerosol pH.
    \param NumberConcentration_aer aerosol concentrations.
  */
  template<class T>
  void Aerosol_SSH<T>::Forward_aer(T current_time,
					  Data<T, 3>& SpecificHumidity_i,
					  Data<T, 3>& Temperature_i,
					  Data<T, 3>& Pressure_i,
					  T delta_t,
					  Data<T, 4>& Concentration,
					  Data<T, 3>&
					  LiquidWaterContent_i,
					  Data<T, 2>& Rain_i,
					  Array<T, 1>& VerticalInterface,
					  Data<T, 5>& Concentration_aer,
					  Data<T, 3>&
					  InCloudWetDepositionFlux,
					  Data<T, 4>&
					  InCloudWetDepositionFlux_aer,
					  Data<T, 3>& pH,
					  Data<T, 4>& NumberConcentration_aer,
					  Data<T, 3>&
					  InCloudWetDepositionFluxNumber_aer)
  {

    int first_index_along_x, last_index_along_x;

    double max_c_aer=0.0;
    Array<T, 1> Index_max(5);
    Index_max=0;
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // This array will store the concentration data with the species
    // axis permuted with the X-axis.
    Array<T, 4> OrdConc;
    Array<T, 5> OrdConc_aer;

    // This array will store the number concentration data
    Array<T, 4> OrdNumConc_aer;
    Array<T, 3> OrdWetDepositionFluxNumber_aer;
	
    Array<T, 3> OrdWetDepositionFlux;
    Array<T, 4> OrdWetDepositionFlux_aer;
    ScatterSlice_x_MPI(Concentration, OrdConc);
    ScatterSlice_x_MPI(Concentration_aer, OrdConc_aer);
    if (option_process_aer["with_number_concentration"])
      {
	ScatterSlice_x_MPI(NumberConcentration_aer, OrdNumConc_aer);
	ScatterSlice_x_MPI(InCloudWetDepositionFluxNumber_aer, OrdWetDepositionFluxNumber_aer);
      }
    ScatterSlice_x_MPI(InCloudWetDepositionFlux, OrdWetDepositionFlux);
    ScatterSlice_x_MPI(InCloudWetDepositionFlux_aer, OrdWetDepositionFlux_aer);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
#else
    first_index_along_x = 0;
    last_index_along_x = Concentration.GetLength(3);
#endif
	
    InCloudWetDepositionFlux.SetZero();
    InCloudWetDepositionFlux_aer.SetZero();
    if (option_process_aer["with_number_concentration"])
      InCloudWetDepositionFluxNumber_aer.SetZero();
    T total_number_in=0.0;
    T total_number_out=0.0;	
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)	\
  firstprivate(first_index_along_x, last_index_along_x)    
#endif
    for (int i = first_index_along_x; i < last_index_along_x; i++)
    {
      for (int k = 0; k < Concentration.GetLength(1); k++)
        for (int j = 0; j < Concentration.GetLength(2); j++)
          {
	    Array<T, 1> Concentration1D(Ns);
            Array<T, 2> Concentration_aer2D(Ns_aer, Nbin_aer);
	    Array<T, 1> NumberConcentration_aer1D(Nbin_aer);
            Array<T, 1> InCloudWetDepositionFlux1D(Ns);
            Array<T, 2> InCloudWetDepositionFlux_aer2D(Ns_aer, Nbin_aer);
	    Array<T, 1> InCloudWetDepositionFluxNumber_aer1D(Nbin_aer);
            Array<T, 1> CurrentVerticalInterface(2);
	    int s, b;

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		InCloudWetDepositionFlux1D(s) = OrdWetDepositionFlux(i, j, s);
#else
		InCloudWetDepositionFlux1D(s) = InCloudWetDepositionFlux(s, j, i);
#endif
	      }

            for (b = 0; b < Nbin_aer ; b++)
              {
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    InCloudWetDepositionFlux_aer2D(s, b) = OrdWetDepositionFlux_aer(i, b, j, s);
#else
                    InCloudWetDepositionFlux_aer2D(s, b) = InCloudWetDepositionFlux_aer(s, b, j, i);
#endif
                  }
              }

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		Concentration1D(s) = OrdConc(i, k, j, s);
#else
		Concentration1D(s) = Concentration(s, k, j, i);
#endif
	      }
            double total_ms=0;
	    double total_nb=0;
            for (b = 0; b < Nbin_aer ; b++)
              {
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    Concentration_aer2D(s, b) = OrdConc_aer(i, b, k, j, s);
#else
                    Concentration_aer2D(s, b) = Concentration_aer(s, b, k, j, i);
		    if(s< Ns_aer-1)
		      total_ms+=Concentration_aer2D(s, b);
#endif
                  }
              }
			
	    if (option_process_aer["with_number_concentration"])
	      {
		for (b = 0; b < Nbin_aer ; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    InCloudWetDepositionFluxNumber_aer1D(b) =
		      OrdWetDepositionFluxNumber_aer(i, j, b);
#else
                    InCloudWetDepositionFluxNumber_aer1D(b) =
		      InCloudWetDepositionFluxNumber_aer(b, j, i);
#endif
                  }
              
		for (b = 0; b < Nbin_aer; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    NumberConcentration_aer1D(b) = OrdNumConc_aer(i, k, j, b);
#else
		    NumberConcentration_aer1D(b) = NumberConcentration_aer(b, k, j, i);
		    total_nb+=NumberConcentration_aer1D(b);
		    total_number_in+= NumberConcentration_aer1D(b);		    
#endif
		  }
	      }
	    else
	      {
		InCloudWetDepositionFluxNumber_aer1D = 0.0;
		NumberConcentration_aer1D = 0.0;
	      }
			
			
            CurrentVerticalInterface(0) = VerticalInterface(k);
            CurrentVerticalInterface(1) = VerticalInterface(k+1);

            Array<T, 1> LiquidWaterContent_1D(Concentration.GetLength(1));
            int nfoglay = 0;
            T lwc_avg = 0.0;
            T heightfog = 0.0;
            int ifog = 0;
            for (int kk = 0; kk < Concentration.GetLength(1); kk++)
              LiquidWaterContent_1D(kk) = LiquidWaterContent_i(kk, j, i)
		* 1000. * Pressure_i(kk, j, i)
		/ 101325.0 * 28.97 / 0.082
		/ Temperature_i(kk, j, i);

            FogSettling(LiquidWaterContent_1D,
                        lwc_cloud_threshold, VerticalInterface,
                        lwc_avg, nfoglay, heightfog);
            if (k < nfoglay) ifog = 1;
            Forward_aer(current_time,
			SpecificHumidity_i(k, j, i), Temperature_i(k, j, i),
			Pressure_i(k, j, i),
			delta_t, Concentration1D,
			LiquidWaterContent_i(k, j, i), Rain_i(j,i),
			CurrentVerticalInterface, Concentration_aer2D,
			InCloudWetDepositionFlux1D,
			InCloudWetDepositionFlux_aer2D,
			pH(k, j, i), lwc_avg, heightfog, ifog,
			NumberConcentration_aer1D,
			InCloudWetDepositionFluxNumber_aer1D);

            
	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		OrdConc(i, k, j, s) = Concentration1D(s);
#else
		Concentration(s, k, j, i) = Concentration1D(s);
#endif
	      }
	    double total_ms_a=0;
	    double total_nb_a=0;
	    for (s = 0; s < Ns_aer; s++)
	      {
                for (b = 0; b < Nbin_aer; b++)
                  {		
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    OrdConc_aer(i, b, k, j, s) = Concentration_aer2D(s, b);
#else
		    Concentration_aer(s, b, k, j, i) = Concentration_aer2D(s, b);

		    if(s< Ns_aer-1)
		    {
		      total_ms_a+=Concentration_aer2D(s, b);
		      if(max_c_aer<Concentration_aer2D(s, b))
		      {
			max_c_aer=Concentration_aer2D(s, b);
			Index_max(0)=s;
			Index_max(1)=b;
			Index_max(2)=k;
			Index_max(3)=j;
			Index_max(4)=i;
		      }
		    }
#endif
                  }
              }

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		OrdWetDepositionFlux(i, j, s) += InCloudWetDepositionFlux1D(s);
#else
		InCloudWetDepositionFlux(s, j, i) += InCloudWetDepositionFlux1D(s);
#endif
	      }

            for (b = 0; b < Nbin_aer ; b++)
              {
                for (s = 0; s < Ns_aer ; s++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    OrdWetDepositionFlux_aer(i, b, j, s) +=
		      InCloudWetDepositionFlux_aer2D(s, b);
#else
                    InCloudWetDepositionFlux_aer(s, b, j, i) +=
		      InCloudWetDepositionFlux_aer2D(s, b);
						
#endif
                  }
              }
			
	    if (option_process_aer["with_number_concentration"])
	      {
		for (b = 0; b < Nbin_aer; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    OrdWetDepositionFluxNumber_aer(i, j, b) +=
		      InCloudWetDepositionFluxNumber_aer1D(b);
#else
		    InCloudWetDepositionFluxNumber_aer(b, j, i) +=
		      InCloudWetDepositionFluxNumber_aer1D(b);
			
#endif
		  }
				
		for (b = 0; b < Nbin_aer ; b++)
		  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
		    OrdNumConc_aer(i, k, j, b) = NumberConcentration_aer1D(b);
#else
		    NumberConcentration_aer(b, k, j, i) = NumberConcentration_aer1D(b);
		    total_nb_a+=NumberConcentration_aer1D(b);
		    total_number_out+= NumberConcentration_aer1D(b);
#endif
		  }
	      }
	  }// loop k,j,i
    }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(OrdConc, Concentration);
    GatherSlice_x_MPI(OrdConc_aer, Concentration_aer);
    if(option_process_aer["with_number_concentration"])
      {
	GatherSlice_x_MPI(OrdNumConc_aer, NumberConcentration_aer);
	GatherSlice_x_MPI(OrdWetDepositionFluxNumber_aer, InCloudWetDepositionFluxNumber_aer);
      }
    GatherSlice_x_MPI(OrdWetDepositionFlux, InCloudWetDepositionFlux);
    GatherSlice_x_MPI(OrdWetDepositionFlux_aer, InCloudWetDepositionFlux_aer);
#endif	
  } // Forward_aer


  template<class T>
  void Aerosol_SSH<T>::Forward(T current_time,
                               T& attenuation,
                               T& humid,
                               T& temperature,
                               T& pressure,
                               Array<T, 1>& source,
                               Array<T, 1>& photolysis,
                               T delta_t,
                               T& attenuation_f,
                               T& humid_f,
                               T& temperature_f,
                               T& pressure_f,
                               Array<T, 1>& source_f,
                               Array<T, 1>& photolysis_f,
                               T& lon, T& lat,
                               Array<T, 1>& concentration,
                               T& liquid_water_content,
                               Array<T, 1>& wet_diameter_aer,
                               Array<T, 2>& concentration_aer,
                               T& ph,
                               Array<T, 1>& number_concentration_aer)
  {

    // Update input data for ssh-aerosol using Polair3d input data.
    UpdateSharedLib(current_time,
                    attenuation,
                    humid,
                    temperature,
                    pressure,
                    delta_t,
                    lon, lat);

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_set_gas_",
                              concentration);


    api.call(_aerosol_so, "api_sshaerosol_updatephoto_");
    

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_set_aero_",
                              concentration_aer);


    if (option_process_aer["with_number_concentration"])
      {
        api.exchange_double_array(_aerosol_so,
                                  "api_sshaerosol_set_aero_num_",
                                  number_concentration_aer);
      }
    else
      {
        // Estimate number concentration from mass
        // Call compute_number
        api.call(_aerosol_so, "api_sshaerosol_compute_number_");
        
      }
    
    // Call gaseous chemistry
    if (with_gas_chemistry)
      api.call(_aerosol_so, "api_sshaerosol_gaschemistry_");


    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_gas_",
                              concentration);

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_aero_",
                              concentration_aer);

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_aero_num_",
                              number_concentration_aer);

  }
  

  template<class T>
  void Aerosol_SSH<T>::Forward_aer(T current_time,
                                   T& specifichumidity,
                                   T& temperature,
                                   T& pressure,
                                   T delta_t,
                                   Array<T, 1>& concentration,
                                   T& liquidwatercontent,
                                   T& rain,
                                   Array<T, 1>& CurrentVerticalInterface,
                                   Array<T, 2>& concentration_aer,
                                   Array<T, 1>& incloudwetdepositionflux,
                                   Array<T, 2>& incloudwetdepositionflux_aer,
                                   T& ph,
                                   T& lwc_avg,
                                   T& heightfog,
                                   int& ifog,
                                   Array<T, 1>& number_concentration_aer,
                                   Array<T, 1>& incloudwetdepositionfluxnumber_aer,
                                   Array<T, 1>& wet_diameter_aer)
  {

    int is_num = option_process_aer["with_number_concentration"] ? 1 : 0;
    int is_fixed_density = option_process_aer["with_fixed_density"] ? 1 : 0;   
    T cloud_water;
    T rain_rate;

    T fixed_density;

    // Density of each grid bin
    Array<T, 1> density_aer_bin;

    // Density of each size section
    Array<T, 1> density_aer_size;
    
    Array<T, 1> bin_diam;

    //! Bins bounds (in µm).
    Array<T, 1> diam_bound;

    //! Averaged particle size in bins (in µm).
    Array<T, 1> size_diam_av;    

    //! Mass bounds (in µm).
    Array<T, 1> mass_bound;

    //! XBF (in µm).
    Array<T, 1> log_bound;

    // 
    Array<int, 2> concentration_index;

    // Index of aerosol speceis
    Array<int, 1> list_species;

    //
    int section_pass;

    // Final time
    T final_time;

    // Set final time.
    final_time = current_time + delta_t;
    
    //! Get aerosol density from SSH-aerosol.
    density_aer_bin.resize(Nbin_aer);
    density_aer_size.resize(Nsize_section_aer);    

    // conversion from kg/m3 to µg/µm3
    fixed_density = FixedDensity_aer * 1.e-9; 
    if (option_process_aer["with_fixed_density"])
      {
        density_aer_bin = fixed_density;
        density_aer_size = fixed_density;
        api.send_double(_aerosol_so,
                        "api_sshaerosol_set_fixed_density_",
                        fixed_density);

      }
    else
      {
        api.call(_aerosol_so,
                 "api_sshaerosol_compute_all_density");
        api.exchange_double_array(_aerosol_so,
                                  "api_sshaerosol_get_density_aer_bin_",
                                  density_aer_bin);
        api.exchange_double_array(_aerosol_so,
                                  "api_sshaerosol_get_density_aer_size_",
                                  density_aer_size);        
      }


    // Aerosol discretization converted m to microm.
    diam_bound.resize(Nsize_section_aer + 1);
    for (int b = 0; b < Nsize_section_aer + 1; b++)
      diam_bound(b) = BinBound_aer(b) * 1.e6;

    mass_bound.resize(Nsize_section_aer + 1);
    log_bound.resize(Nsize_section_aer + 1);    
    for (int b = 0; b < Nsize_section_aer + 1; b++)
      {
        if (b == Nsize_section_aer)
          mass_bound(b) = density_aer_size(b - 1)* cst_pi6 * pow(diam_bound(b), 3);
        else
          mass_bound(b) = density_aer_size(b) * cst_pi6 * pow(diam_bound(b), 3);

        if (mass_bound(b) <= 0.0)
          throw ("Error: zero or negative mass_bound");
        log_bound(b) = log10(mass_bound(b));
      }
    
    // Get averaged particle size in each bin (in µm).
    size_diam_av.resize(Nsize_section_aer);
    for (int b = 0; b < Nsize_section_aer; b++)
      size_diam_av(b) = sqrt(diam_bound(b) * diam_bound(b + 1));
    
    //! Cloud liquid water content (g/m^3)
    cloud_water = liquidwatercontent * 1000.0
      * pressure / 101325.0 * 28.97 / Pr / temperature;


    //! Aqueous chemistry is taken into account
    //! if liquid water content is larger than a minimum value.
    if (cloud_water >= lwc_cloud_threshold)
      {
        if (option_process_aer["with_in_cloud_scavenging"])
          rain_rate = rain;
        else
          rain_rate = 0.0;


        concentration_index.resize(2, Nbin_aer);
        api.exchange_int_array(_aerosol_so,
                               "api_sshaerosol_get_concentration_index_",
                               concentration_index);

        // Get list_species from ssh-aerosol
        list_species.resize(Ns_aer);
        api.exchange_int_array(_aerosol_so,
                               "api_sshaerosol_get_list_species_",
                               list_species);        

        //
        section_pass = api.recv_int(_aerosol_so,
                                    "api_sshaerosol_get_section_pass_");
        

        // Call VSRM aqueous chemistry module.
        if (aqueous_module == "vsrm")
          {
            // Estimate number concentration from mass
            if (!option_process_aer["with_number_concentration"])
              {
                // Call compute_number
                api.call(_aerosol_so, "api_sshaerosol_compute_number_");

                // Get number concentration from ssh-aerosol
                api.exchange_double_array(_aerosol_so,
                                          "api_sshaerosol_get_aero_num_",
                                          number_concentration_aer);
              }

            // Call VSRM
            _vsrmchem(&Ns, &Ns_aer, &Nbin_aer, &Nsize_section_aer, &Ns_cloud_interact,
                      density_aer_bin.data(), &fixed_density, BinBound_aer.data(),
                      size_diam_av.data(), log_bound.data(),
                      concentration.data(), concentration_aer.data(),
                      &specifichumidity, &pressure, &temperature,
                      &cloud_water, &current_time, &final_time, &rain_rate, &ph,
                      incloudwetdepositionflux.data(),
                      incloudwetdepositionflux_aer.data(),
                      incloudwetdepositionfluxnumber_aer.data(),
                      &is_num,
                      number_concentration_aer.data(), &conserving_mass_tolerance,
                      &species_index_cloud_interact.front(),
                      concentration_index.data(), list_species.data(),
                      &Ns_inorganic_aer, &Ns_organic_aer, &section_pass, &iredist,
                      &is_fixed_density);
          }
        // Call a simplified aqueous chemistry module.        
        else if (aqueous_module == "simple")
          {
	    if (Ns_cloud_interact != 15)
	      {
		throw string( "be careful, the gas_species_cloud_interact list in the species file must, for now, contain exactly 15 elements. The list must respect this order : NH3 HNO3 HCL SO2 H2O2 HCHO HONO O3 OH HO2 NO3 NO NO2 PAN H2SO4" );
	      }
            // Estimate number concentration from mass
            if (!option_process_aer["with_number_concentration"])
              {
                // Call compute_number
                api.call(_aerosol_so, "api_sshaerosol_compute_number_");

                // Get number concentration from ssh-aerosol
                api.exchange_double_array(_aerosol_so,
                                          "api_sshaerosol_get_aero_num_",
                                          number_concentration_aer);
              }
            
            _simple_aqueous_module
              (&Ns, &Ns_aer, &Nbin_aer, &Nsize_section_aer,
               &Ns_cloud_interact,
               density_aer_bin.data(), &fixed_density, BinBound_aer.data(),
               size_diam_av.data(), log_bound.data(),
               concentration.data(), concentration_aer.data(),
               &specifichumidity, &pressure, &temperature,
               &cloud_water, &current_time, &final_time, &rain_rate, &ph,
               incloudwetdepositionflux.data(),
               incloudwetdepositionflux_aer.data(),
               incloudwetdepositionfluxnumber_aer.data(),
               &is_num,
               number_concentration_aer.data(), &conserving_mass_tolerance,
               &species_index_cloud_interact.front(),
               concentration_index.data(), list_species.data(),
               &Ns_inorganic_aer, &Ns_organic_aer, &section_pass, &iredist,
               &is_fixed_density);
          }
        // Call ssh-aerosol
        else
          {

            
            // Send concentration arrays to ssh-aerosol.
            api.exchange_double_array(_aerosol_so,
                                      "api_sshaerosol_set_gas_",
                                      concentration);

            api.exchange_double_array(_aerosol_so,
                                      "api_sshaerosol_set_aero_",
                                      concentration_aer);

            if (option_process_aer["with_number_concentration"])
              {
                api.exchange_double_array(_aerosol_so,
                                          "api_sshaerosol_set_aero_num_",
                                          number_concentration_aer);
              }
            else
              {
                // Estimate number concentration from mass
                // Call compute_number
                api.call(_aerosol_so, "api_sshaerosol_compute_number_");
                
              }
            
            
            // Call aerosols dynamic
            api.call(_aerosol_so, "api_sshaerosol_aerodyn_");

          }
      }
    // Low cloud_water
    else
      {
        // Send concentration arrays to ssh-aerosol.
        api.exchange_double_array(_aerosol_so,
                                  "api_sshaerosol_set_gas_",
                                  concentration);

        api.exchange_double_array(_aerosol_so,
                                  "api_sshaerosol_set_aero_",
                                  concentration_aer);

        if (option_process_aer["with_number_concentration"])
          {
            api.exchange_double_array(_aerosol_so,
                                      "api_sshaerosol_set_aero_num_",
                                      number_concentration_aer);
          }
        else
          {
            // Estimate number concentration from mass
            // Call compute_number
            api.call(_aerosol_so, "api_sshaerosol_compute_number_");
          }
        
        // Call aerosols dynamic
        api.call(_aerosol_so, "api_sshaerosol_aerodyn_");
      }
   

    // Return new gas-phase concentrations.
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_gas_",
                              concentration);

    // Return new aerosol mass concentrations.
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_aero_",
                              concentration_aer);

    // Return new aerosol number concentrations.
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_aero_num_",
                              number_concentration_aer);


    // Return wet diameter.
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_wet_diameter_",
                              wet_diameter_aer);


    // Output from SSH-aerosol.
    // Please use compiling option ssh-output=yes
    // if you run 0-D test cases. 
#ifdef WRITE_SSH_OUTPUT            
        api.call(_aerosol_so, "api_sshaerosol_output_");
#endif

  }

  //! Checks whether a given field is required by this module.
  /*! Checks whether a given field must be available in the underlying model
    for this module to perform properly.
    \param field the field name.
    \return True is the field is required, false otherwise.
  */
  template<class T>
  bool Aerosol_SSH<T>::IsRequired(string field)
  {
    if (field == "wet_diameter_aer") return true;
    if (field == "rain") return true;
    return false;
  }


  //! Checks whether a given field is computed by this module.
  /*! Checks whether a given field is computed by this module and updated in
    the underlying model.
    \param field the field name.
    \return True is the field is computed, false otherwise.
  */
  template<class T>
  bool Aerosol_SSH<T>::IsComputed(string field)
  {
    if (field == "pH") return true;
    return false;
  }

  template<class T>
  void Aerosol_SSH<T>::FogSettling(Array<T, 1>& LiquidWaterContent,
                                   T& lwc_cloud_threshold,
                                   Array<T, 1>& VerticalInterface,
                                   T& lwc_avg,
                                   int& nfoglay,
                                   T& heightfog)
  {
    // Fog settling
    int size_lwc = LiquidWaterContent.shape()[0];
	
    if (LiquidWaterContent(0) > lwc_cloud_threshold)
      {
	int indok = 1;
	for (int k = 0; k < size_lwc; k++)
	  {
	    if ((LiquidWaterContent(k) > lwc_cloud_threshold) && (indok == 1))
	      {
				
		heightfog = VerticalInterface(k+1);
		nfoglay = k+1;
		lwc_avg += LiquidWaterContent(k);
			
	      }
	    else
	      indok = 0;
	  }
		
	if (nfoglay > 0) 
	  lwc_avg /= nfoglay;
      }
   
  }
 

  //! Update input data for ssh-aerosol shared library.
  template<class T>
  void Aerosol_SSH<T>::UpdateSharedLib(T current_time,
                                       T attenuation,
                                       T humidity,
                                       T temperature,
                                       T pressure,
                                       T delta_t,
                                       T lon, T lat)
  {
    
    T rh;
    _compute_relative_humidity(&humidity, &temperature,
                               &pressure, &rh);
    // cout << "++ Computed Relative humidity: " << rh << endl;
    
    api.send_double(_aerosol_so,
                    "api_sshaerosol_set_current_t_",
                    current_time);

    api.send_double(_aerosol_so,
                    "api_sshaerosol_set_attenuation_",
                    attenuation);
    
    api.send_double(_aerosol_so,
                    "api_sshaerosol_set_humidity_",
                    humidity);

    api.send_double(_aerosol_so,
                    "api_sshaerosol_set_relhumidity_",
                    rh);
    
    api.send_double(_aerosol_so,
                    "api_sshaerosol_set_temperature_",
                    temperature);

    api.send_double(_aerosol_so,
                    "api_sshaerosol_set_pressure_",
                    pressure);

    api.send_double(_aerosol_so,
                    "api_sshaerosol_set_dt_",
                    delta_t);

    T lonlat[2] = {lon, lat};
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_set_lonlat_",
                              lonlat);

  }



  //! Get Ns_aer from ssh-aerosol.
  template<class T>
  Array<T, 1>
  Aerosol_SSH<T>::GetMassDensityLayers()
  {
    return ssh_mass_density_layers;
  }
  
  //! Get Ns_aer_nolayer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNsAerNoLayer()
  {
    return Ns_aer_nolayer;
  }

  
  //! Get Ns_aer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNsAer()
  {
    return Ns_aer;
  }


  //! Get N_layer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNLayer()
  {
    return N_layer;
  }
  
  //! Get i_hydrophilic from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetIHydrophilic()
  {
    return i_hydrophilic;
  }
  
  //! Get Nbin_aer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNbinAer()
  {
    return Nbin_aer;
  }

  //! Get nsize_section_aer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNSizeSectionAer()
  {
    return Nsize_section_aer;
  }


  
  //! Whether externally mixed.
  template<class T>
  bool Aerosol_SSH<T>::IsExternalComposition()
  {
    return with_external_composition;
  }  

  //! Get index_groups from ssh-aerosol.
  template<class T>
  Array<int, 1> Aerosol_SSH<T>::GetIndexGroups()
  {
    return index_groups;
  }

  //! Get index_groups_ext from ssh-aerosol.
  template<class T>
  Array<int, 1> Aerosol_SSH<T>::GetIndexGroupsExt()
  {
    return index_groups_ext;
  }  
  

  //! Get Ngroup_aer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNgroup_aer()
  {
    return Ngroup_aer;
  }

  //! Get Ncomposition_aer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNcomposition_aer()
  {
    return Ncomposition_aer;
  }

  //! Get Nfraction_aer from ssh-aerosol.
  template<class T>
  int Aerosol_SSH<T>::GetNfraction_aer()
  {
    return Nfraction_aer;
  }  

  //! Get index_groups_ext from ssh-aerosol.
  template<class T>
  Array<T, 3> Aerosol_SSH<T>::GetCompositionBounds()
  {

    discretization_composition_conv.resize(Ncomposition_aer, Ngroup_aer, 2);
    
    // Array conversion from (2, Ngroup_aer, Ncompostion_aer)
    // to (Ncomposition_aer, Ngroup_aer, 2)
    for (int i = 0; i < Ncomposition_aer; i++)
      for (int j = 0; j < Ngroup_aer; j++)
        for (int k = 0; k < 2; k++)
          discretization_composition_conv(i, j, k) = discretization_composition(k, j, i);
    
    return discretization_composition_conv;
  }  

  //! Get aerosol_type from ssh-aerosol.
  template<class T>
  Array<int, 1> Aerosol_SSH<T>::GetAerosolType()
  {
    return aerosol_type;
  }


  //! Get aerosol species name from SSH-aerosol.
  template<class T>
  Array<string, 1> Aerosol_SSH<T>::GetAerosolSpecName()
  {
    return aerosol_spec_name;
  }
  
  
//   //! Performs an integration over one time step.
//   /*!
//     \param current_time starting time in seconds.
//     \param Attenuation_i cloud attenuation coefficients at the beginning of
//     the time step.
//     \param SpecificHumidity_i specific humidity at the beginning of the time
//     step.
//     \param Temperature_i temperature at the beginning of the time step.
//     \param Pressure_i pressure at the beginning of the time step.
//     \param Source_i volume sources at the beginning of the time step.
//     \param PhotolysisRate_i photolysis rates at the beginning of the time
//     step.
//     \param next_time time at the end of the time step.
//     \param Attenuation_f cloud attenuation coefficients at the end of the
//     time step.
//     \param SpecificHumidity_f specific humidity at the end of the time step.
//     \param Temperature_f temperature at the end of the time step.
//     \param Pressure_f pressure at the end of the time step.
//     \param Source_f volume sources at the end of the time step.
//     \param PhotolysisRate_f photolysis rates at the end of the time step.
//     \param Longitude longitudes.
//     \param Latitude latitudes.
//     \param Concentration gas concentrations.
//     \param LiquidWaterContent_i air liquid water content at the beginning of
//     the time step.
//     \param WetDiameter_aer Aerosol wet diameters at the beginning of the time
//     step.
//     \param Concentration_aer aerosol concentrations.
//   */
  template<class T>
  template<class ClassModel>
  void Aerosol_SSH<T>
  ::InitDistribution(ClassModel& Model)
    
  {

    InitDistribution(Model.GetConcentration(),
                     Model.GetConcentration_aer(),
                     Model.GetNumberConcentration_aer());

  }

  //! Initialize the aerosol distribution for MUNICH.
  /*!

   */
  template<class T>
  void Aerosol_SSH<T>
  ::InitDistribution(Data<T, 2>& Concentration,
                     Data<T, 3>& Concentration_aer,
                     Data<T, 2>& NumberConcentration_aer)
  {
    
    int first_index_along_x, last_index_along_x;
    first_index_along_x = 0;
    last_index_along_x = Concentration.GetLength(1);

    for (int i = first_index_along_x; i < last_index_along_x; i++)
      {
        Array<T, 1> Concentration1D(Ns);
        Array<T, 2> Concentration_aer2D(Ns_aer_nolayer, Nsize_section_aer);            
        Array<T, 2> InitConcentration_aer(Ns_aer, Nbin_aer);

        Array<T, 1> NumberConcentration_aer1D(Nsize_section_aer);            
        NumberConcentration_aer1D = 0.0;

        Array<T, 1> InitNumberConcentration_aer(Nbin_aer);
            
        int s, b;

        for (s = 0; s < Ns; s++)
          {
            Concentration1D(s) = Concentration(s, i);
          }

        int Ncomposition = Nbin_aer / Nsize_section_aer;
        int isize = 0;
        int ibin = 0;
        for (b = 0; b < Nbin_aer ; b++)              
          {
            int isp = 0;
            int ilay = 0;
            // Take the first bin of each size
            if (ibin == 0)
              {
                for (s = 0; s < Ns_aer ; s++)
                  {
                    InitConcentration_aer(s, b) = Concentration_aer(s, b, i);
                        
                    // Inorganic species and water
                    if (s < (Ns_inorganic_aer + Ns_inert_aer) or (s == Ns_aer - 1))
                      {
                        Concentration_aer2D(isp, isize) = InitConcentration_aer(s, b);
                        ++isp;
                      }
                    // Organic species in the layers
                    else
                      {
                        // Take concentrations of the first layer.
                        if (ilay == 0)
                          {
                            Concentration_aer2D(isp, isize) = InitConcentration_aer(s, b);
                            ++isp;
                          }
                        ++ilay;

                        if (ilay == (N_layer + i_hydrophilic))
                          ilay = 0;
                      }
                  }
                ++isize;
              }
            ++ibin;
            if (ibin == Ncomposition)
              ibin = 0;
          }
            
            
        if (option_process_aer["with_number_concentration"])
          {
            ibin = 0;
            isize = 0;
            for (b = 0; b < Nbin_aer; b++)
              {
                // Take the first bin of each size
                if (ibin == 0)
                  {
                    InitNumberConcentration_aer(b) =
                      NumberConcentration_aer(b, i);
                        
                    NumberConcentration_aer1D(isize) = InitNumberConcentration_aer(b);
                    ++isize;
                  }
                ++ibin;
                if (ibin == Ncomposition)
                  ibin = 0;
              }
          }
            

        InitDistribution(Concentration1D,
                         Concentration_aer2D,
                         NumberConcentration_aer1D,
                         InitConcentration_aer,
                         InitNumberConcentration_aer);
            

        // Return concentrations from init_distribution.
        for (s = 0; s < Ns; s++)
          {
            Concentration(s, i) = Concentration1D(s);
          }

        for (s = 0; s < Ns_aer; s++)
          {
            for (b = 0; b < Nbin_aer; b++)
              {
                Concentration_aer(s, b, i) = InitConcentration_aer(s, b);
              }
          }
            

        if (option_process_aer["with_number_concentration"])
          for (b = 0; b < Nbin_aer ; b++)
            {
              NumberConcentration_aer(b, i) = InitNumberConcentration_aer(b);
            } 

        
      }
    
  } // InitDistribution



  
  template<class T>
  void Aerosol_SSH<T>
  ::InitDistribution(Data<T, 4>& Concentration,
                     Data<T, 5>& Concentration_aer,
                     Data<T, 4>& NumberConcentration_aer)
  {
    
    int first_index_along_x, last_index_along_x;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // This array will store the concentration data with the species
    // axis permuted with the X-axis.
    Array<T,4> OrdConc;
    Array<T,5> OrdConc_aer;
    Array<T,4> OrdNumConc_aer;
    ScatterSlice_x_MPI(Concentration, OrdConc);
    ScatterSlice_x_MPI(Concentration_aer, OrdConc_aer);
    if (option_process_aer["with_number_concentration"])
      ScatterSlice_x_MPI(NumberConcentration_aer, OrdNumConc_aer);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
#else
    first_index_along_x = 0;
    last_index_along_x = Concentration.GetLength(3);
#endif
    
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)	\
  firstprivate(first_index_along_x, last_index_along_x)
#endif
    for (int i = first_index_along_x; i < last_index_along_x; i++)
      for (int k = 0; k < Concentration.GetLength(1); k++)
        for (int j = 0; j < Concentration.GetLength(2); j++)
	  {
            Array<T, 1> Concentration1D(Ns);
            Array<T, 2> Concentration_aer2D(Ns_aer_nolayer, Nsize_section_aer);
            Array<T, 2> InitConcentration_aer(Ns_aer, Nbin_aer);
	    Array<T, 1> NumberConcentration_aer1D(Nsize_section_aer);
            NumberConcentration_aer1D = 0.0;

	    Array<T, 1> InitNumberConcentration_aer(Nbin_aer);
            
            int s, b;

	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		Concentration1D(s) = OrdConc(i, k, j, s);
#else
		Concentration1D(s) = Concentration(s, k, j, i);
#endif
	      }

            int Ncomposition = Nbin_aer / Nsize_section_aer;
            int isize = 0;
            int ibin = 0;
            for (b = 0; b < Nbin_aer ; b++)              
              {
                int isp = 0;
                int ilay = 0;
                // Take the first bin of each size
                if (ibin == 0)
                  {
                    for (s = 0; s < Ns_aer ; s++)
                      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
                        // Rank of indices must match the one defined in
                        // ScatterSlice_x_MPI.
                        InitConcentration_aer(s, b) = OrdConc_aer(i, b, k, j, s);
#else
                        InitConcentration_aer(s, b) = Concentration_aer(s, b, k, j, i);
#endif
                        // Inorganic species and water
                        if (s < (Ns_inorganic_aer + Ns_inert_aer) or (s == Ns_aer - 1))
                          {
                            Concentration_aer2D(isp, isize) = InitConcentration_aer(s, b);
                            ++isp;
                          }
                        // Organic species in the layers
                        else
                          {
                            // Take concentrations of the first layer.
                            if (ilay == 0)
                              {
                                Concentration_aer2D(isp, isize) = InitConcentration_aer(s, b);
                                ++isp;
                              }
                            ++ilay;

                            if (ilay == (N_layer + i_hydrophilic))
                              ilay = 0;
                          }
                      }
                    ++isize;
                  }
                ++ibin;
                if (ibin == Ncomposition)
                  ibin = 0;
              }
            
            
	    if (option_process_aer["with_number_concentration"])
              {
                ibin = 0;
                isize = 0;
                for (b = 0; b < Nbin_aer; b++)
                  {
                    // Take the first bin of each size
                    if (ibin == 0)
                      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
                        // Rank of indices must match the one defined in
                        // ScatterSlice_x_MPI.

                        if (isnan(OrdNumConc_aer(i, k, j, b)))
                          cout <<"isnan "<< b <<" "<< i <<" "<<  j<<" " <<  k;

                        InitNumberConcentration_aer(b) =
                          OrdNumConc_aer(i, k, j, b);
#else
                        InitNumberConcentration_aer(b) =
                          NumberConcentration_aer(b, k, j, i);
#endif
                        
                        NumberConcentration_aer1D(isize) = InitNumberConcentration_aer(b);
                        ++isize;
                      }
                    ++ibin;
                    if (ibin == Ncomposition)
                      ibin = 0;
                  }
              }
            

            InitDistribution(Concentration1D,
                             Concentration_aer2D,
                             NumberConcentration_aer1D,
                             InitConcentration_aer,
                             InitNumberConcentration_aer);

            // Return concentrations from init_distribution.
	    for (s = 0; s < Ns; s++)
	      {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		// Rank of indices must match the one defined in
		// ScatterSlice_x_MPI.
		OrdConc(i, k, j, s) = Concentration1D(s);
#else
		Concentration(s, k, j, i) = Concentration1D(s);
#endif
	      }

	    for (s = 0; s < Ns_aer; s++)
              {
                for (b = 0; b < Nbin_aer; b++)
                  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		    // Rank of indices must match the one defined in
		    // ScatterSlice_x_MPI.
                    OrdConc_aer(i, b, k, j, s) = InitConcentration_aer(s, b);
#else
                    Concentration_aer(s, b, k, j, i) = InitConcentration_aer(s, b);
#endif
                  }
              }

	    if (option_process_aer["with_number_concentration"])
	      for (b = 0; b < Nbin_aer ; b++)
		{
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
		  // Rank of indices must match the one defined in
		  // ScatterSlice_x_MPI.
		  OrdNumConc_aer(i, k, j, b) = InitNumberConcentration_aer(b);
#else
                  NumberConcentration_aer(b, k, j, i) = InitNumberConcentration_aer(b);
#endif
		} 


          }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(OrdConc, Concentration);
    GatherSlice_x_MPI(OrdConc_aer, Concentration_aer);
    if(option_process_aer["with_number_concentration"])
      GatherSlice_x_MPI(OrdNumConc_aer, NumberConcentration_aer);
#endif
    
  } // InitDistribution
  

  template<class T>
  void Aerosol_SSH<T>
  ::InitDistribution(Array<T, 1>& Concentration1D,
                     Array<T, 2>& Concentration_aer2D,
                     Array<T, 1>& NumberConcentration_aer1D,
                     Array<T, 2>& InitConcentration_aer,
                     Array<T, 1>& InitNumberConcentration_aer)
  {

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_set_gas_",
                              Concentration1D);

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_set_init_bin_mass_",
                              Concentration_aer2D);

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_set_init_bin_number_",
                              NumberConcentration_aer1D);
    api.call(_aerosol_so,
             "api_sshaerosol_init_distributions_");

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_gas_",
                              Concentration1D);
            
    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_aero_",
                              InitConcentration_aer);

    api.exchange_double_array(_aerosol_so,
                              "api_sshaerosol_get_aero_num_",
                              InitNumberConcentration_aer);


    // Output from SSH-aerosol.
    // Please use compiling option ssh-output=yes
    // if you run 0-D test cases. 
#ifdef WRITE_SSH_OUTPUT            
    /* InitOutput */
    api.call(_aerosol_so, "api_sshaerosol_initoutput_");
                
    /* Report */
    api.call(_aerosol_so, "api_sshaerosol_report_");
        
    /* Output */
    api.call(_aerosol_so, "api_sshaerosol_output_");
#endif
    
  }
  
} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SSH_CXX
#endif


