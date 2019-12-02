#ifndef POLYPHEMUS_FILE_MODELS_STREETNETWORKCHEMISTRY_CXX

//////////////
// INCLUDES //

#include "StreetNetworkChemistry.hxx"

// INCLUDES //
//////////////
  
namespace Polyphemus
{

  template<class T, class ClassChemistry>
  const T StreetNetworkChemistry<T, ClassChemistry>::pi = acos(-1);

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds the model. Nothing else is performed.
  */
  template<class T, class ClassChemistry>
  StreetNetworkChemistry<T, ClassChemistry>::StreetNetworkChemistry():
    StreetNetworkTransport<T>()
  {
  }


  //! Main constructor.
  /*!
    \param config_file configuration filename.
  */
  template<class T, class ClassChemistry>
  StreetNetworkChemistry<T, ClassChemistry>::StreetNetworkChemistry(string config_file):
    StreetNetworkTransport<T>(config_file)
  {
  }
  

  //! Destructor.
  template<class T, class ClassChemistry>
  StreetNetworkChemistry<T, ClassChemistry>::~StreetNetworkChemistry()
  {
    this->ClearStreetVector();
    this->ClearIntersectionVector();
  }

  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    species lists, options (especially which processes are included) and the
    paths to data input-files.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::ReadConfiguration()
  {
    StreetNetworkTransport<T>::ReadConfiguration();

    /*** Options ***/

    this->config.SetSection("[options]");

    this->config.PeekValue("With_chemistry",
			   this->option_process["with_chemistry"]);

    this->config.PeekValue("Option_chemistry", option_chemistry);

    this->config.PeekValue("With_photolysis",
			   this->option_process["with_photolysis"]);

    CheckConfiguration();

  }

  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::Init()
  {
    //    StreetNetworkTransport<T>::Init();
    this->SetCurrentDate(this->Date_min);
    InitStreet();
    this->InitIntersection();

    Allocate();

    this->ComputeStreetAngle();

    if (this->option_process["with_chemistry"])
      InitChemistry();
   
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Initialize the number of sources to parallelize.
    BaseModuleParallel::InitSource(this->GetNStreet());

    //
    BaseModuleParallel::BuildPartition_source();
#endif
    
  }

  //! Streets initialization.
  /*!
   */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::InitStreet()
  {
    Array<T, 1> init_conc(this->Ns);
    init_conc = 0.0;
    for (int i = 0; i < this->total_nstreet; ++i)
      {
        Street<T>* street = 
          new Street<T>(this->id_street(i), this->begin_inter(i), this->end_inter(i),
                        this->length(i), this->width(i), this->height(i),
                        this->Ns, Nr_photolysis);

        this->StreetVector.push_back(street);
      }
    this->current_street = this->StreetVector.begin();
  }

  //! Checks the configuration.
  /*! 
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::CheckConfiguration()
  {
    // The configuration-file path is the field "Data_description" in the main
    // configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");
    // Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    // Photolysis rates.
    if (this->option_process["with_photolysis"])
      this->input_files["photolysis_rates"].ReadFiles(data_description_file,
						      "photolysis_rates");
    else
      this->input_files["photolysis_rates"].Empty();
    for (map<string, string>::iterator i
	   = this->input_files["photolysis_rates"].Begin();
	 i != this->input_files["photolysis_rates"].End(); i++)
      photolysis_reaction_list.push_back(i->first);
    Nr_photolysis = int(photolysis_reaction_list.size());

    // Reads data description.
    if (this->option_process["with_photolysis"])
      {
	data_description_stream.SetSection("[photolysis_rates]");
	photolysis_date_min = data_description_stream.PeekValue("Date_min");
	data_description_stream.PeekValue("Delta_t", "> 0",
					  photolysis_delta_t);
	data_description_stream.PeekValue("Ndays", "> 0",
					  Nphotolysis_days);
	data_description_stream.PeekValue("Time_angle_min",
					  photolysis_time_angle_min);
	data_description_stream.PeekValue("Delta_time_angle", "> 0",
					  photolysis_delta_time_angle);
	data_description_stream.PeekValue("Ntime_angle", "> 0",
					  Nphotolysis_time_angle);
	data_description_stream.PeekValue("Latitude_min",
					  photolysis_latitude_min);
	data_description_stream.PeekValue("Delta_latitude", "> 0",
					  photolysis_delta_latitude);
	data_description_stream.PeekValue("Nlatitude", "> 0",
					  Nphotolysis_latitude);
	data_description_stream.Find("Altitudes");
	split(data_description_stream.GetLine(), altitudes_photolysis);
      }
  }

  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::Allocate()
  {
    StreetNetworkTransport<T>::Allocate();

    /*** Photolysis rates ***/
    
    if (this->option_process["with_photolysis"])
      {
	GridR_photolysis = RegularGrid<T>(Nr_photolysis);
    
	Grid_time_angle_photolysis =
	  RegularGrid<T>(photolysis_time_angle_min,
			 photolysis_delta_time_angle,
			 Nphotolysis_time_angle);
	Grid_latitude_photolysis = RegularGrid<T>(photolysis_latitude_min,
						  photolysis_delta_latitude,
						  Nphotolysis_latitude);

	Nphotolysis_z = int(altitudes_photolysis.size());
	GridZ_photolysis = RegularGrid<T>(Nphotolysis_z);
	for (unsigned int i = 0; i < altitudes_photolysis.size(); i++)
	  GridZ_photolysis(i) = to_num<T>(altitudes_photolysis[i]);

      }
  }


  //! Initializes photolysis rates.
  /*! The rates are computed on the basis of raw photolysis-rates read in
    files. The raw photolysis-rates depend upon the day, the time angle, the
    latitude and the altitude.
    \param date the date at which photolysis rates are needed.
    \param Rates (output) the photolysis rates.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::InitPhotolysis(Date date)
  {
    Array<T, 1> photolysis_rate(Nr_photolysis);

    int i, j, k;
    T time_angle;
    int angle_in, j_in, k_in;
    T alpha_angle, alpha_y, alpha_z;
    T one_alpha_angle, one_alpha_y, one_alpha_z;
    int nb_days;

    // Relevant step in input files.
    int day = int(T(date.GetSecondsFrom(photolysis_date_min))
		  / 86400. / photolysis_delta_t + .5);

    if (day >= Nphotolysis_days)
      throw string("There are not enough available data for photolysis ")
	+ string("rates. Missing days.");

    Data<T, 3> FileRates(Grid_time_angle_photolysis,
			 Grid_latitude_photolysis, GridZ_photolysis);

    // Interpolation.
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        T altitude = 2.0;
        T latitude = street->GetLatitude();
        T longitude = street->GetLongitude();

        // Loop over photolysis reactions.
        for (int r = 0; r < Nr_photolysis; r++)
          {
            string filename =
              this->input_files["photolysis_rates"](photolysis_reaction_list[r]);

            FormatBinary<float> format;
            format.ReadRecord(filename, day, FileRates);
          
            // Along z.
            k_in = 0;
            while (k_in < Nphotolysis_z - 1 && GridZ_photolysis(k_in)
                   < altitude)
              k_in++;
            if (k_in == Nphotolysis_z - 1
                && GridZ_photolysis(k_in) < altitude)
              throw string("There are not enough available data for ")
                + string("photolysis rates. Missing levels.");
            if (k_in > 0)
              k_in--;
            alpha_z = (altitude -
                       GridZ_photolysis(k_in))
              / (GridZ_photolysis(k_in + 1) - GridZ_photolysis(k_in));

            // Along y (latitude).
            j_in = int((latitude
                        - photolysis_latitude_min)
                       / photolysis_delta_latitude);
            alpha_y = (latitude
                       - photolysis_latitude_min
                       - T(j_in) * photolysis_delta_latitude)
              / photolysis_delta_latitude;

            // Time angle.
            time_angle = T(date.GetNumberOfSeconds()) / 3600.
              - 12. + longitude / 15.;
            nb_days = int(time_angle / 24.);
            time_angle = abs(time_angle - 24. * T(nb_days));
            if (time_angle > 12.)
              time_angle = 24. - time_angle;
		
            angle_in = int((time_angle - photolysis_time_angle_min)
                           / photolysis_delta_time_angle);
            alpha_angle = (time_angle - photolysis_time_angle_min
                           - T(angle_in) * photolysis_delta_time_angle)
              / photolysis_delta_time_angle;

            one_alpha_angle = 1. - alpha_angle;
            one_alpha_y = 1. - alpha_y;
            one_alpha_z = 1. - alpha_z;

            if (angle_in >= Nphotolysis_time_angle - 1)
              photolysis_rate(r) = 0.;
            else
              photolysis_rate(r) =
                one_alpha_z * one_alpha_y * one_alpha_angle
                * FileRates(angle_in, j_in, k_in)
                + alpha_z * one_alpha_y * one_alpha_angle
                * FileRates(angle_in, j_in, k_in + 1)
                + one_alpha_z * alpha_y * one_alpha_angle
                * FileRates(angle_in, j_in + 1, k_in)
                + alpha_z * alpha_y * one_alpha_angle
                * FileRates(angle_in, j_in + 1, k_in + 1)
                + one_alpha_z * one_alpha_y * alpha_angle
                * FileRates(angle_in + 1, j_in, k_in)
                + alpha_z * one_alpha_y * alpha_angle
                * FileRates(angle_in + 1, j_in, k_in + 1)
                + one_alpha_z * alpha_y * alpha_angle
                * FileRates(angle_in + 1, j_in + 1, k_in)
                + alpha_z * alpha_y * alpha_angle
                * FileRates(angle_in + 1, j_in + 1, k_in + 1);
          }
        street->SetPhotolysisRate(photolysis_rate);
      }
    
  }



  //! Performs one step forward.
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::Forward()
  {
    this->Transport();
    
    if (this->option_process["with_chemistry"])
      {
	if (this->option_process["with_stationary_hypothesis"])
	    Chemistry();
	else
	    ComputeStreetConcentrationNoStationary(); // without stationarity
      }

    this->SetStreetConcentration();

    this->AddTime(this->Delta_t);
    this->step++;

  }

  //LL Remove stationary regime
  //! Compute the concentrations in the street-canyon using the flux balance equation.
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::ComputeStreetConcentrationNoStationary()
  {
    // MPI implementation.
    // 'For' loop on the street segments is parallelized.
    int first_index_along_source, last_index_along_source;
    
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    
    int nstreet = this->StreetVector.size();

    // The following arrays need to be used only for MPI.
    Array<T, 2> concentration_mpi;
    Array<T, 2> massflux_roof_to_background_mpi;
    Array<T, 2> street_quantity_delta_mpi;
  
    concentration_mpi.resize(nstreet, this->Ns);
    massflux_roof_to_background_mpi.resize(nstreet, this->Ns);
    street_quantity_delta_mpi.resize(nstreet, this->Ns);
    
    concentration_mpi = 0.0;
    massflux_roof_to_background_mpi = 0.0;
    street_quantity_delta_mpi = 0.0;
    
    // Scatter buffers for Street-type arrays to all processes.
    BaseModuleParallel::ScatterSlice_source_MPI(concentration_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(massflux_roof_to_background_mpi);
    BaseModuleParallel::ScatterSlice_source_MPI(street_quantity_delta_mpi);

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
	Array<T, 1> concentration_array(this->Ns);
	Array<T, 1> concentration_array_tmp(this->Ns);
	Array<T, 1> init_concentration_array(this->Ns);
	Array<T, 1> background_concentration_array(this->Ns);
	Array<T, 1> new_concentration_array(this->Ns);
	Array<T, 1> emission_rate_array(this->Ns);
	Array<T, 1> inflow_rate_array(this->Ns);
	Array<T, 1> deposition_flux_array(this->Ns);

	concentration_array = 0.0;
	concentration_array_tmp = 0.0;
	init_concentration_array = 0.0;
	background_concentration_array = 0.0;
	new_concentration_array = 0.0;
	emission_rate_array = 0.0;
	inflow_rate_array = 0.0;
	deposition_flux_array = 0.0;

        T transfer_velocity = street->GetTransferVelocity(); // m/s
	T temp = transfer_velocity * street->GetWidth() * street->GetLength(); // m3/s
        T outgoing_flux = street->GetOutgoingFlux(); // m3/s
        T inflow_flux = street->GetIncomingFlux();
        T street_volume = street->GetHeight() *
          street->GetWidth() * street->GetLength(); // m3
            
	for (int s = 0; s < this->Ns; ++s)
          {
	    concentration_array(s) = street->GetStreetConcentration(s);
	    init_concentration_array(s) = street->GetStreetConcentration(s);
	    background_concentration_array(s) = street->GetBackgroundConcentration(s);
	    emission_rate_array(s) = street->GetEmission(s);
	    inflow_rate_array(s) = street->GetInflowRate(s);
	    deposition_flux_array(s) = street->GetDepositionRate() * street_volume; // to update
	  }

        T sub_delta_t_init, sub_delta_t;
	Date current_date_tmp = this->current_date;

        // Initialize a time-step of solver. 
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
		   deposition_flux_array);
        
	Date next_date = this->current_date;
	Date next_date_tmp = this->current_date;
	next_date.AddSeconds(this->Delta_t);
	next_date_tmp.AddSeconds(sub_delta_t_init);

	while (current_date_tmp < next_date)
	  {
            // YK : removed because it is already read.
	    // //! Get street concentrations.
	    // for (int s = 0; s < this->Ns; ++s)
	    //   concentration_array(s) = street->GetStreetConcentration(s);

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
                                   deposition_flux_array,
                                   new_concentration_array,
                                   sub_delta_t_init,
                                   street->GetStreetID());

	      }
	    else if (this->option_method == "Rosenbrock")
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
					  deposition_flux_array,
					  new_concentration_array,
					  sub_delta_t_init,
					  inflow_flux);
	      }
	    else 
	      throw string("Error: numerical method not chosen.");

            
            //! sub_delta_t_init corresponds to the time step being incremented
            //! sub_delta_t corresponds to the time step that may be used in the next iteration
	    //! Calculates the new sub_delta_t for the next iteration
	    T sub_delta_t_max = next_date.GetSecondsFrom(next_date_tmp);
	    StreetNetworkTransport<T>::
	      AdaptTimeStep(new_concentration_array,
			    concentration_array_tmp,
			    sub_delta_t_init,
			    this->sub_delta_t_min,
			    sub_delta_t_max,
			    sub_delta_t);

	    //! Chemical reactions
	    if (this->option_process["with_chemistry"])
	      {
		Array<T, 1> photolysis_rate(Nr_photolysis);

		T attenuation_ = street->GetAttenuation();
		T specific_humidity_ = street->GetSpecificHumidity();
		T pressure_ = street->GetPressure();
		T temperature_ = street->GetTemperature();
		T longitude_ = street->GetLongitude();
		T latitude_ = street->GetLatitude();

		for (int r = 0; r < Nr_photolysis; r++)
		  photolysis_rate(r) = street->GetPhotolysisRate(r);

		Chemistry(current_date_tmp,
			  new_concentration_array,
			  attenuation_,
			  specific_humidity_,
			  pressure_,
			  temperature_,
			  longitude_,
			  latitude_,
			  photolysis_rate,
			  sub_delta_t_init);
	      }
	    
	    //! Set the new concentrations.
	    for (int s = 0; s < this->Ns; ++s)
              street->SetStreetConcentration(new_concentration_array(s), s);

	    //! Actualises current_time_tmp
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
      }

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI

    // Gather data from processes to Process 0
    // and cast all data from Process 0 to all other processes.
    BaseModuleParallel::GatherSlice_source_MPI(concentration_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(massflux_roof_to_background_mpi);
    BaseModuleParallel::GatherSlice_source_MPI(street_quantity_delta_mpi);

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
      }
#endif
    
  }


  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::InitData()
  {
    StreetNetworkTransport<T>::InitData();

    /*** Meteo for the chemistry on the streets ***/
    string filename;
    if (this->option_process["with_local_data"])
      {
        attenuation.resize(this->Nt_meteo, this->total_nstreet);
        filename = this->input_files["meteo"]("Attenuation");
        InitData(filename, attenuation);

        specific_humidity.resize(this->Nt_meteo, this->total_nstreet);
        filename = this->input_files["meteo"]("SpecificHumidity");
        InitData(filename, specific_humidity);

        pressure.resize(this->Nt_meteo, this->total_nstreet);
        filename = this->input_files["meteo"]("SurfacePressure");
        InitData(filename, pressure);

        temperature.resize(this->Nt_meteo, this->total_nstreet);
        filename = this->input_files["meteo"]("SurfaceTemperature");
        InitData(filename, temperature);
      }

  }

  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::InitData(string input_file, 
                                                           Array<T, 2>& input_data)
  {
    StreetNetworkTransport<T>::InitData(input_file, input_data);
  }

  //! Method called at each time step to initialize the model.
  /*!
    \note Empty method.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::InitStep()
  {
    StreetNetworkTransport<T>::InitStep();
    
    InitPhotolysis(this->current_date);

    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        if (this->option_process["with_local_data"])
          {
            //! Set the meteo data
            street->SetMeteoChemistry(attenuation(this->meteo_index, st),
                                      specific_humidity(this->meteo_index, st),
                                      pressure(this->meteo_index, st),
                                      temperature(this->meteo_index, st));
          }
	st += 1;
      }

  }


  //! Chemistry.
  /*!
    \param 
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>::Chemistry()
  {

    Array<T, 1> source(this->Ns);
    Array<T, 1> photolysis_rate(Nr_photolysis);
    Array<T, 1> concentration_array(this->Ns);

    source = 0.0;
    
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin(); 
	 iter != this->StreetVector.end(); iter++)
      {
	Street<T>* street = *iter;

	//! Get the concentrations.
	for (int s = 0; s < this->Ns; ++s)
	  concentration_array(s) = street->GetStreetConcentration(s);

	T attenuation_ = street->GetAttenuation();
	T specific_humidity_ = street->GetSpecificHumidity();
	T pressure_ = street->GetPressure();
	T temperature_ = street->GetTemperature();
	T longitude_ = street->GetLongitude();
	T latitude_ = street->GetLatitude();

	for (int r = 0; r < Nr_photolysis; r++)
	  photolysis_rate(r) = street->GetPhotolysisRate(r);

#ifndef POLYPHEMUS_WITH_AEROSOL_MODULE
	Chemistry_.Forward(T(this->current_date.GetNumberOfSeconds()),
			     attenuation_, specific_humidity_,
			     temperature_, pressure_, source,
			     photolysis_rate,
			     T(this->next_date.
			       GetSecondsFrom(this->current_date)),
			     attenuation_, specific_humidity_,
			     temperature_, pressure_, source,
			     photolysis_rate, longitude_,
			     latitude_, concentration_array);

	for (int s = 0; s < this->Ns; ++s)
	  street->SetStreetConcentration(concentration_array(s), s);
#endif	
      }
  }

  //! Chemistry in no stationary regime.
  /*!
    \param 
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>
  ::Chemistry(Date current_date_tmp,
	      Array<T, 1>& concentration_array,
	      T attenuation_,
	      T specific_humidity_,
	      T pressure_,
	      T temperature_,
	      T longitude_,
	      T latitude_,
	      Array<T, 1> photolysis_rate,
	      T sub_delta_t)
  {
    Array<T, 1> source(this->Ns);
    source = 0.0;
    Chemistry_.Forward(T(current_date_tmp.GetNumberOfSeconds()),
		       attenuation_, specific_humidity_,
		       temperature_, pressure_, source,
		       photolysis_rate,
		       sub_delta_t,
		       attenuation_, specific_humidity_,
		       temperature_, pressure_, source,
		       photolysis_rate, longitude_,
		       latitude_, concentration_array);
  }

  /*! Initializes background concentrations and photolysis parameters.
    \param meteo ConfigStream instance through which background concentrations
    and photolysis rates may be read.
  */
  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>
  ::InitChemistry()
  {
    InitPhotolysis(this->current_date);
    
    Chemistry_.Init(*this);
  }


  template<class T, class ClassChemistry>
  void StreetNetworkChemistry<T, ClassChemistry>
  ::SetStreetAdditionalMeteo(int street_index, T attenuation,
                             T specific_humidity, T pressure,
                             T temperature)
  {
    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = this->StreetVector.begin();
         iter != this->StreetVector.end(); iter++)
      {
        if (st == street_index)
          {
            Street<T>* street = *iter;
            street->SetMeteoChemistry(attenuation, specific_humidity, 
                              pressure, temperature);
          }
        ++st;
      }

  }

  /*!
  */
  template<class T, class ClassChemistry>
  bool StreetNetworkChemistry<T, ClassChemistry>
  ::WithChemistry()
  {
    return this->option_process["with_chemistry"];
  }

  /*!
  */
  template<class T, class ClassChemistry>
  string StreetNetworkChemistry<T, ClassChemistry>
  ::GetChemicalMechanism()
  {
    return option_chemistry;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKCHEMISTRY_CXX
#endif
