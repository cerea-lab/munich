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

    this->config.SetSection("[street]");
    this->config.PeekValue("Transfert_parameterization",
                           "Sirane | Schulte", option_transfer);
    this->config.PeekValue("Mean_wind_speed_parameterization",
                           "Sirane | Lemonsu", option_ustreet);
    this->config.PeekValue("Numerical_method_parameterization",
			   "ETR | Rosenbrock", option_method);
    this->config.PeekValue("Building_height_wind_speed_parameterization",
                           "Sirane | Macdonald", option_uH);
    this->config.PeekValue("With_horizontal_fluctuation",
			   this->option_process["with_horizontal_fluctuation"]);

    this->config.PeekValue("Minimum_Street_Wind_Speed", ">= 0.1", ustreet_min);

    this->config.PeekValue("With_local_data",
			   this->option_process["with_local_data"]);
    this->config.PeekValue("With_stationary_hypothesis",
			   this->option_process["with_stationary_hypothesis"]);

    this->config.PeekValue("Sub_delta_t_min", sub_delta_t_min);

    
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

    if (this->option_process["with_local_data"])
      CheckConfiguration();

    ReadStreetData();

    Compute_z0_d_city();

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif
    if (rank == 0)
      DisplayConfiguration();

  }

  //! Display the configuration
  template<class T>
  void StreetNetworkTransport<T>::DisplayConfiguration()
  {
    cout << "Transfert_parameterization: " << option_transfer << endl;
    cout << "Mean_wind_speed_parameterization: " << option_ustreet << endl;
    cout << "Numerical_method_parameterization: " << option_method << endl;
    cout << "Building_height_wind_speed_parameterization: " << option_uH << endl;
  }
  
  //! Allocates memory.
  /*! Allocates the grids and the concentration Data.
   */
  template<class T>
  void StreetNetworkTransport<T>::Allocate()
  {
    GridST2D = RegularGrid<T>(total_nstreet);
    GridS2D = RegularGrid<T>(this->Ns);

    /*** Street concentration ***/
    StreetConcentration.Resize(GridS2D, GridST2D);
    StreetConcentration.SetZero();
  }

  //! Model initialization.
  /*! It reads the configuration, allocates memory and reads the values of the
    fields at the beginning of the simulation.
  */
  template<class T>
  void StreetNetworkTransport<T>::Init()
  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif
    
    this->SetCurrentDate(this->Date_min);
    InitStreet();
    InitIntersection();

    Allocate();

    ComputeStreetAngle();
  }


  //! Street source initialization.
  template<class T>
  void StreetNetworkTransport<T>::InitData(string input_file, Array<T, 2>& input_data)
  {
    if (is_num(input_file))
      {
	T value = to_num<T>(input_file);
        if (value < 0.0)
          throw string("Negative value is not allowed. ") + input_file;
	input_data = value;
	return;
      }
    FormatBinary<float>().Read(input_file, input_data);    

  }

  //! Input data initialization.
  template<class T>
  void StreetNetworkTransport<T>::InitData()
  {
    string filename;

    /*** Meteo for the streets ***/
    if (this->option_process["with_local_data"])
      {
        wind_direction_arr.resize(Nt_meteo, total_nstreet);
        wind_speed_arr.resize(Nt_meteo, total_nstreet);
        pblh_arr.resize(Nt_meteo, total_nstreet);
        ust_arr.resize(Nt_meteo, total_nstreet);
        lmo_arr.resize(Nt_meteo, total_nstreet);

        filename = this->input_files["meteo"]("WindDirection");
        InitData(filename, wind_direction_arr);
        filename = this->input_files["meteo"]("WindSpeed");
        InitData(filename, wind_speed_arr);
        filename = this->input_files["meteo"]("PBLH");
        InitData(filename, pblh_arr);
        filename = this->input_files["meteo"]("UST");
        InitData(filename, ust_arr);
        filename = this->input_files["meteo"]("LMO");
        InitData(filename, lmo_arr);
    
        /*** Meteo for the intersections ***/
        wind_direction_inter.resize(Nt_meteo, nintersection);
        wind_speed_inter.resize(Nt_meteo, nintersection);
        pblh_inter.resize(Nt_meteo, nintersection);
        ust_inter.resize(Nt_meteo, nintersection);
        lmo_inter.resize(Nt_meteo, nintersection);

        filename = this->input_files["meteo"]("WindDirectionInter");
        InitData(filename, wind_direction_inter);
        filename = this->input_files["meteo"]("WindSpeedInter");
        InitData(filename, wind_speed_inter);
        filename = this->input_files["meteo"]("PBLHInter");
        InitData(filename, pblh_inter);
        filename = this->input_files["meteo"]("USTInter");
        InitData(filename, ust_inter);
        filename = this->input_files["meteo"]("LMOInter");
        InitData(filename, lmo_inter);
      }

    /*** Emission data ***/

    emission.resize(Ns_emis, Nt_emis, total_nstreet);
    for (int s = 0; s < Ns_emis; s++)
      {
        filename = this->input_files["emission"](species_list_emis[s]);
        Array<T, 2> emission_tmp(Nt_emis, total_nstreet);

        InitData(filename, emission_tmp);
        for (int t = 0; t < Nt_emis; t++)
          for (int st = 0; st < total_nstreet; st++)
            emission(s, t, st) = emission_tmp(t, st);

      }

    /*** Background concentration ***/

    if (this->option_process["with_local_data"])
      {
        background_concentration.resize(Ns_background, Nt_background, total_nstreet);
        for (int s = 0; s < Ns_background; s++)
          {
            filename = this->input_files["background_concentration"](species_list_background[s]);
            Array<T, 2> background_tmp(Nt_background, total_nstreet);
            InitData(filename, background_tmp);
            for (int t = 0; t < Nt_background; t++)
              for (int st = 0; st < total_nstreet; st++)
                background_concentration(s, t, st) = background_tmp(t, st);
          }
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
    for (int i = 0; i < total_nstreet; ++i)
      {
        StreetStream.GetLine(line);
        v = split(line, ";");
	id_street(i) = to_num<T>(v[0]);
        begin_inter(i) = to_num<T>(v[1]);
        end_inter(i) = to_num<T>(v[2]);
        length(i) = to_num<T>(v[3]);
        width(i) = to_num<T>(v[4]);
        height(i) = to_num<T>(v[5]);
	
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
    
  }

  //! Computation of the displacement and roughness heights of the Macdonald wind profile
  template<class T>
  void StreetNetworkTransport<T>::Compute_z0_d_city()
  {
    cout << "Mean street length: " << Mean_length << "m" << endl;
    cout << "Mean street width: " << Mean_width << "m" << endl;
    cout << "Mean buildings height: " << Mean_height << "m" << endl;
    T building_width = Mean_width; // hyp: W street = W buildings
    cout << "Mean buildings width: " << building_width << endl;
    T Cd_building = 1.2; // Drag coefficient for buildings
    T betta = 0.55; // Correction coefficient betta = 1.0 for staggered arrays and betta = 0.55 for square array (Macdonald et al, 1998)
    T A = 3.59; // Empiric coefficient A = 4.43 for staggered arrays and A = 3.59 for the square arrays (Macdonald et al, 1998)   

    T Af = Mean_length * Mean_height; // frontal area of obstacles
    T Ap = Mean_length * building_width; // plan area of obstacles
    T At = Mean_length * (Mean_width + building_width); // lot area of obstacles

    T Lambdaf = Af / At; // frontal area density of obstacles
    T Lambdap = Ap / At; // plan area density of obstacles

    d_city = Mean_height * (1.0 + pow(A,-Lambdap) * (Lambdap - 1.0)); // displacement height of city
    T temp_dH = 1.0 - (d_city / Mean_height);
    z0_city = Mean_height * (temp_dH * exp(-pow(0.5 * betta * Cd_building / pow(karman, 2.0) * temp_dH * Lambdaf,-0.5))); //roughness length of city
    cout << "Hauteur de deplacement d_city : " << d_city << "m" << endl;
    cout << "Hauteur de rugosite z0_city : " << z0_city << "m" << endl;

  }

  //! Streets initialization.
  /*!
   */
  template<class T>
  void StreetNetworkTransport<T>::InitStreet()
  {

    Array<T, 1> init_conc(this->Ns);
    init_conc = 0.0;
    int nr_photolysis = 24;
    for (int i = 0; i < total_nstreet; ++i)
      {
        Street<T>* street = 
          new Street<T>(id_street(i), begin_inter(i), end_inter(i),
                        length(i), width(i), height(i),
                        this->Ns, nr_photolysis);
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
    Array<T, 1> emission_rate(this->Ns);
    Array<T, 2> street_emission;
    Date date_min;
    T time_distance, delta_t;

    int background_index = 9999;

    date_min = this->input_files["emission"].GetDateMin();
    time_distance = this->current_date.GetSecondsFrom(date_min);
    delta_t = this->input_files["emission"].GetDelta_t();
    int emis_index = 9999;
    if (time_distance < 0.0)
      throw string("Error: missing emission input data: ") 
        + "the current date is " + to_str(this->current_date) 
        + "\n but the beginning date for the street emission data is " + to_str(date_min);
    else
      {
        emis_index = time_distance / delta_t;
        if (emis_index >= emission.shape()[1])
          throw("Error: missing emission input data.");
      }

    if (this->option_process["with_local_data"])
      {    
        date_min = this->input_files["meteo"].GetDateMin();
        time_distance = this->current_date.GetSecondsFrom(date_min);
        delta_t = this->input_files["meteo"].GetDelta_t();
        meteo_index = time_distance / delta_t;
        if (meteo_index >= wind_direction_arr.shape()[0])
          throw("Error: missing meteo input data.");

        date_min = this->input_files["background_concentration"].GetDateMin();
        time_distance = this->current_date.GetSecondsFrom(date_min);
        delta_t = this->input_files["background_concentration"].GetDelta_t();
        background_index = time_distance / delta_t;
        if (background_index >= background_concentration.shape()[1])
          throw("Error: missing background concentration data.");
      }

    int st = 0;
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin();
         iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;

        //! Set the emission data
        emission_rate = 0.0; 
        for (int i = 0; i < this->Ns; ++i)
          for (int j = 0; j < Ns_emis; ++j)
            if (this->species_list[i] == species_list_emis[j])
	      { 
		emission_rate(i) = emission(j, emis_index, st);
	      }
	//throw;
        street->SetEmission(emission_rate);

        if (this->option_process["with_local_data"])
          {
            //! Set the meteo data
            street->SetMeteo(wind_direction_arr(meteo_index, st),
                             wind_speed_arr(meteo_index, st),
                             pblh_arr(meteo_index, st),
                             ust_arr(meteo_index, st),
                             lmo_arr(meteo_index, st));

            //! Set the background concentration data
            for (int i = 0; i < this->Ns; ++i)
              for (int j = 0; j < Ns_background; ++j)
                if (this->species_list[i] == species_list_background[j])
                  street->SetBackgroundConcentration(background_concentration(j, background_index, st), i);
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
            intersection->SetMeteo(wind_direction_inter(emis_index, inter),
                                   wind_speed_inter(emis_index, inter),
                                   pblh_inter(emis_index, inter),
                                   ust_inter(emis_index, inter),
                                   lmo_inter(emis_index, inter));
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
    /*** Meteo ***/

    //! Check the meteo data for the streets.
    if (this->input_files["meteo"]("WindDirection").empty())
      throw "WindDirection is needed but no input data file was provided.";
    if (this->input_files["meteo"]("WindSpeed").empty())
      throw "WindSpeed is needed but no input data file was provided.";
    if (this->input_files["meteo"]("PBLH").empty())
      throw "PBLH is needed but no input data file was provided.";
    if (this->input_files["meteo"]("UST").empty())
      throw "UST is needed but no input data file was provided.";

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


  //! Performs one step forward.
  template<class T>
  void StreetNetworkTransport<T>::Forward()
  {
    Transport();

    SetStreetConcentration();

    this->AddTime(this->Delta_t);
    this->step++;
  }

  //! Calls the functions for the transport.
  template<class T>
  void StreetNetworkTransport<T>::Transport()
  {

    is_stationary = false;

    //! Compute the wind speed in the street-canyon.
   
    if (option_uH == "Macdonald")
   	Compute_Macdonald_uH();

    ComputeUstreet();

    ComputeSigmaW();

    ComputeTransferVelocity();

    ComputeWindDirectionFluctuation();

    SetInitialStreetConcentration();

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

	//! Compute the concentrations in the street-canyon.
	if (this->option_process["with_chemistry"] == false)
	  ComputeStreetConcentrationNoStationary();
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
        T street_volume = street->GetVolume(); // m3        
        
        bool is_stationary_local = true;
        for (int s = 0; s < this->Ns; ++s)
          {
            T street_conc = street->GetStreetConcentration(s); // ug/m3
            T init_street_conc = street->GetInitialStreetConcentration(s); // ug/m3
	    T emission_rate = street->GetEmission(s); //  ug/s
            T inflow_rate = street->GetInflowRate(s); // ug/s
            T deposition_rate = street->GetDepositionRate(); // 1/s
            T deposition_flux = deposition_rate * street_volume; // m3/s
            T conc_bg = street->GetBackgroundConcentration(s); // ug/m3

            //! Compute the new concentrations.
            T street_conc_new = 0.0;
            if ((temp + outgoing_flux + deposition_flux) == 0.0)
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
                (temp + outgoing_flux + deposition_flux);

            //! Set the minimum possible concentrations.
            street_conc_new = max(street_conc_new, 0.0);

            //! Check the stationarity
            T delta_concentration = abs(street_conc_new - street_conc);
            if ((street_conc != 0.0) and 
                (delta_concentration > delta_concentration_min))
              is_stationary_local = false;

            //! Set the new concentrations.
            street->SetStreetConcentration(street_conc_new, s);

            T massflux_roof_to_background;
            if (is_stationary_local)
              {
                massflux_roof_to_background = temp * (street_conc_new - conc_bg); // ug/s

                street->SetMassfluxRoofToBackground(massflux_roof_to_background, s);

                T conc_delta = street_conc_new - init_street_conc;
                T street_quantity_delta = conc_delta * street_volume; // ug
                street->SetStreetQuantityDelta(street_quantity_delta, s);
              }

          }

        //! True if stationarity is obtained for a street.
        street->SetStationary(is_stationary_local);
      }
  }

  //LL: Remove stationary regime
  //! Compute the concentrations in the street-canyon using the flux balance equation.
  template<class T>
  void StreetNetworkTransport<T>::ComputeStreetConcentrationNoStationary()
  {
    
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
        T street_volume = street->GetVolume();
	T inflow_flux = street->GetIncomingFlux();

	for (int s = 0; s < this->Ns; ++s)
          {
	    concentration_array(s) = street->GetStreetConcentration(s);
	    init_concentration_array(s) = street->GetStreetConcentration(s);
	    background_concentration_array(s) = street->GetBackgroundConcentration(s);
	    emission_rate_array(s) = street->GetEmission(s);
	    inflow_rate_array(s) = street->GetInflowRate(s);
	    deposition_flux_array(s) = street->GetDepositionRate() * street_volume;
	  }

	T sub_delta_t_init, sub_delta_t;

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
		 deposition_flux_array);

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
				 deposition_flux_array,
				 new_concentration_array,
				 sub_delta_t_init,
                                 street->GetStreetID());
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
					deposition_flux_array,
					new_concentration_array,
					sub_delta_t_init,
					inflow_flux);
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
        T sigma_w = street->GetSigmaW();
        if (option_transfer == "Sirane")
          street->SetTransferVelocity(sigma_w / (sqrt(2.0) * pi));
        else if (option_transfer == "Schulte")
          {
            T aspect_ratio = street->GetHeight() / street-> GetWidth();
            const T beta = 0.45;
            T velocity = beta * sigma_w * (1.0 / (1.0 + aspect_ratio));
            street->SetTransferVelocity(velocity);
          }               
      }
  }
  
  //! Compute uH with the method of Macdonald et al (1998)
  template<class T>
  void StreetNetworkTransport<T>::Compute_Macdonald_uH()
  {
    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
    {
        Street<T>* street = *iter;
        T h = street->GetHeight();
        T u_zref = street->GetWindSpeed(); // wind speed at the reference altitude (m/s)
        T zref = h + 17; // reference altitude (m), to change according to our input data
        ustar_macd = u_zref * karman / log((zref - d_city)/z0_city);
        street->SetStreetUstar(ustar_macd);
        if (h < (d_city + z0_city))
        {
          uH_macd = 0.0;
        }
        else
        {
          uH_macd = ustar_macd / karman * log((h - d_city)/z0_city);
        }
    }

    for (typename vector<Intersection<T>* >::iterator iter = IntersectionVector.begin(); iter != IntersectionVector.end(); iter++)
    {
    	Intersection<T>* intersection = *iter;
	intersection->SetIntersectionUST(ustar_macd);
    }
  }


  //! Compute the wind speed in the street-canyon.
  //! Option "Sirane" based on Soulhac et al. (2008)
  //! Option "Lemonsu" based on Lemonsu et al. (2004), and Cherin et al. (2015)
  template<class T>
  void StreetNetworkTransport<T>::ComputeUstreet()
  {
    T phi, delta_i;
    //! Wall roughness length from WRF/UCM
    T z0_build = 0.001;
    T beta;
    Array<T, 1> z, ustreet_z;
    T ustreet;
    int nz = 10;
    z.resize(nz);
    ustreet_z.resize(nz);

    //LL: print data
    bool isIsolated;
    if(total_nstreet > 30)
      isIsolated = false;
    else
      isIsolated = true;

    for (typename vector<Street<T>* >::iterator iter = StreetVector.begin(); iter != StreetVector.end(); iter++)
      {
        Street<T>* street = *iter;
        T h = street->GetHeight();
        T w = street->GetWidth();
        T ang = street->GetStreetAngle();
        T ustar = street->GetStreetUstar();
        T u_zref = street->GetWindSpeed(); // wind speed at the reference altitude (m/s)
        T wind_direction = street->GetWindDirection();

        delta_i = min(h, w / 2.0);
        phi = abs(wind_direction - ang);
        
        T solutionC = ComputeBesselC(z0_build, delta_i);

        if (option_ustreet == "Lemonsu")
          {
            if (( h == 0.0) or (w == 0.0))
              throw string("Street height or width is zero. ") + to_str(h) + " " + to_str(w);
            beta = h / (2.0 * w);
            ustreet = 0.0;
            if (option_uH == "Sirane")
              {
                u_h = ComputeUH(solutionC, ustar);
              }
            else if (option_uH == "Macdonald")
              {
                u_h = uH_macd;
              }
            else
              throw("Wrong option given. Choose Macdonald or Sirane");
            for (int k = 0; k < nz; ++k)
              {
                z(k) = h / (nz - 1) * k;
                ustreet_z(k) = u_h * abs(cos(phi)) * exp(beta * (z(k) / h - 1.0));
                ustreet += ustreet_z(k); 
              }
            ustreet /= nz;
          }
        else if (option_ustreet == "Sirane")
          {
            T u_h_sirane = ComputeUH(solutionC, ustar);
            T alpha = log(delta_i / z0_build);
            T beta = exp(solutionC / sqrt(2.0) * (1.0 - h / delta_i));
            T temp1 = pow(delta_i, 2.0) / (h * w);
            T temp2 = 2.0 * sqrt(2.0) / solutionC * (1.0 - beta);
            T temp3 = 1.0 - pow(solutionC, 2.0) / 3.0 + pow(solutionC, 4.0) / 45.0;
            T temp4 = beta * (2.0 * alpha - 3.0) / alpha;
            T temp5 = (w / delta_i - 2.0) * (alpha - 1.0) / alpha;
            ustreet = u_h_sirane * abs(cos(phi)) * temp1 * 
              (temp2 * temp3 + temp4 + temp5);
          }
        else
          throw("Wrong option given. Choose Lemonsu or Sirane");

        street->SetStreetWindSpeed(max(ustreet, ustreet_min));
      }
  }

  //! Compute U_H in the formulation for the wind speed in the street-canyon.
  /*!
    \param c solution of the Bessel function.
    \param ustar Friction velocity.
  */
  template<class T>
  T StreetNetworkTransport<T>::ComputeUH(T c, T ustar)
  {
    T j1C = j1(c);
    T y0C = y0(c);
    T y1C = y1(c);
    T j0C = j0(c);
    T term1 = pi / (sqrt(2.0) * pow(karman, 2.0) * c);
    T term2 = y0C - (j0C * y1C) / j1C;
    return ustar * sqrt(term1 * term2);
  }
  
  //! Compute the solution C of the equation (1) in Soulhac et al. (2011).
  //! y1(C) and j1(C): the Bessel function of the first kind of order 0 and 1 (see 
  //! http://www.gnu.org/software/libc/manual/html_node/Special-Functions.html)
  /*!
    \param z0_build  Wall roughness length
    \param delta_i
  */
  template<class T>
  T StreetNetworkTransport<T>::ComputeBesselC(T z0_build, T delta_i)
  {
    T gamma = 0.577;
    int nc = 100;
    T maxC = 2.0;
    T step = maxC / nc;
    Array<T, 1> temp(nc);
    Array<T, 1> listC(nc);
    T tempC = 0.0;
    T solutionC;
    for (int i = 0; i < nc; ++i)
      {
        tempC += step;
        listC(i) = tempC;
        T y1C = y1(tempC);
        T j1C = j1(tempC);
        temp(i) = abs(2.0 / tempC * exp(pi / 2.0 * y1C / j1C - gamma) - 
                      z0_build / delta_i);
      }

    int index;
    T minValue;
    GetMin(temp, nc, minValue, index);
    solutionC = listC(index);
    T threshold = 0.001;
    if (minValue > threshold)
      throw string("Fail to find a solution. Please adjust maxC (current value: ")
        + to_str(maxC) + ").";
    return solutionC;
  }

  //! Returns the minimum value of array
  template<class T>
  void StreetNetworkTransport<T>::GetMin(Array<T, 1> arr, int length, T& minimum, int& index)
  {
    minimum = arr(0);
    index = 0;
    for (int i = 1; i < length; ++i)
      if (minimum > arr(i))
        {
          minimum = arr(i);
          index = i; 
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
    
    output_stream.PeekValue("Output_dir", output_dir);

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

  //LL: Remove stationary hypothesis
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
		     Array<T, 1>& new_concentration_array,
		     const T sub_delta_t,
                     const int street_id)
  {
    T dtetr;
    Array<T, 1> dc1dt(this->Ns);
    Array<T, 1> dc2dt(this->Ns);
    dc1dt = 0.0;
    dc2dt = 0.0;

    
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
	       dc1dt);

    
    
    //Update concentrations
    for (int s = 0; s < this->Ns; s++)
      {
	concentration_array_tmp(s) = concentration_array(s) + dc1dt(s)*sub_delta_t;
	concentration_array_tmp(s) = max(concentration_array_tmp(s), 0.0);
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
	       dc2dt);

    
    dtetr = sub_delta_t/2.0;
    for (int s = 0; s < this->Ns; s++)
      {
	new_concentration_array(s) = concentration_array(s) + (dc1dt(s) + dc2dt(s))*dtetr;
	new_concentration_array(s) = max(new_concentration_array(s), 0.0);
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
	     const Array<T, 1> deposition_flux_array)

  {
    Array<T, 1> dcdt(this->Ns);
    dcdt = 0.0;
    Array<T, 1> sub_delta_t_sp(this->Ns);
    sub_delta_t_sp = 0.0;
    Array<T, 1> concentration_array_tmp(this->Ns);
    concentration_array_tmp = 0.0;
    
    if(this->current_date == this->Date_min)
      {
	Array<T, 1> dcdt0(this->Ns);
	dcdt0 = 0.0;
	Array<T, 1> dcdtmin(this->Ns);
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
		   dcdt0);
	
	for (int s = 0; s < this->Ns; s++)
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
		   dcdtmin);
	for (int s = 0; s < this->Ns; s++)
	  sub_delta_t_sp(s) = dcdt0(s)/dcdtmin(s) * sub_delta_t_min;
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
		   dcdt);
	for (int s = 0; s < this->Ns; s++)
	  {
	    T tmp = concentration_array(s) * dcdt(s);
	    if (tmp != 0.0)
	      {
		T tscale = concentration_array(s)/abs(dcdt(s));
		sub_delta_t_sp(s) = min(tscale, this->Delta_t);
		sub_delta_t_sp(s) = max(tscale, sub_delta_t_min);
	      }
	  }
      }
    sub_delta_t = min(sub_delta_t_sp);
    sub_delta_t = max(sub_delta_t, sub_delta_t_min);
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
	       Array<T, 1>& dcdt)
  {
    for (int s = 0; s < this->Ns; s++)
      dcdt(s) = (emission_rate_array(s) + inflow_rate_array(s) - outgoing_flux*concentration_array(s) - deposition_flux_array(s)*concentration_array(s) - temp*(concentration_array(s) - background_concentration_array(s)))/street_volume;
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
			    const T inflow_flux)
  {
    Array<T, 1> dc1dt(this->Ns); // Differential term dc/dt at initial time (sub time step)
    Array<T, 1> dc2dt(this->Ns); // Differential term dc/dt at final time (sub time step)
    Array<T, 1> k1(this->Ns);
    Array<T, 1> k2(this->Ns);
    T gamma;
    dc1dt = 0.0;
    dc2dt = 0.0;
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
	       dc1dt);
    // Compute the Jacobian at initial time
    Array<T, 1> J(this->Ns);
    Jacobian(inflow_flux,
	     outgoing_flux,
	     temp,
	     deposition_flux_array,
	     J,
	     street_volume);
  
    // Compute k1
    for (int s = 0; s < this->Ns; s++)
      k1(s) = dc1dt(s) / (1 - gamma * sub_delta_t * J(s));

    // Compute temporary concentrations at final time
    for (int s = 0; s < this->Ns; s++)
      {
	concentration_array_tmp(s) = concentration_array(s) + k1(s) * sub_delta_t;
	concentration_array_tmp(s) = max(concentration_array_tmp(s), 0.0);
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
	       dc2dt);
    
    // Compute k2
    for (int s = 0; s < this->Ns; s++)
      k2(s) = (dc2dt(s) - 2 * k1(s)) / (1 - gamma * sub_delta_t * J(s));
    
    // Compute concentrations at final time
    for (int s = 0; s < this->Ns; s++)
      {
	new_concentration_array(s) = concentration_array(s) + (3 * k1(s) + k2(s)) * sub_delta_t / 2;
	new_concentration_array(s) = max(new_concentration_array(s), 0.0);
      }
  }
  
  // Compute Jacobian matrix (1x1) by Euler
  template<class T>
  void StreetNetworkTransport<T>
  ::Jacobian(const T inflow_flux,
	     const T outgoing_flux,
	     const T temp,
	     const Array<T, 1> deposition_flux_array,
	     Array<T, 1>& J,
	     const T street_volume)
  {
    for (int s = 0; s < this->Ns; s++)
      {
	J(s) = (inflow_flux - temp - outgoing_flux - deposition_flux_array(s)) / street_volume;
	//J(s) = 0.01;
	//cout<<J(s)<<endl;
      }
  }
  

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETNETWORKTRANSPORT_CXX
#endif
