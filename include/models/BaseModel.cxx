// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
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


#ifndef POLYPHEMUS_FILE_MODELS_BASEMODEL_CXX


#include "BaseModel.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*! Builds the model. Nothing else is performed: the configuration is not
    read and nothing is allocated.
    \param config_file configuration file.
  */
  template<class T>
  BaseModel<T>::BaseModel(string config_file):
    file_config(config_file), config(config_file), current_time(0.), step(0),
    backward(false)
  {
    field_species["Concentration"] = &species_list;
    field_species["Concentration_aer"] = &species_list_aer;

    // May be useful in case of external initializations.
    Ns = 0;
    Ns_aer = 0;
    Nbin_aer = 0;
  }


  //! Destructor.
  template<class T>
  BaseModel<T>::~BaseModel()
  {
  }


  ///////////////////
  // CONFIGURATION //
  ///////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    the species list and display options.
  */
  template<class T>
  void BaseModel<T>::ReadConfiguration()
  {
    /*** Dates ***/

    config.SetSection("[domain]");

    Date_min = config.PeekValue("Date_min");
    config.PeekValue("Delta_t", "> 0", Delta_t);
    config.PeekValue("Nt", "positive", Nt);

    /*** Species ***/

    config.PeekValue("Species", file_species);
    // Opens the file that describes species.
    ConfigStream species_stream(file_species);
    // Section "[species]" contains all species names.
    species_stream.SetSection("[species]");
    while (!species_stream.IsEmpty())
      species_list.push_back(species_stream.GetElement());
    Ns = int(species_list.size());

    /*** Spatial extent ***/

    config.PeekValue("Vertical_levels", file_vertical_levels);
    config.PeekValue("Nz", "positive", Nz);
    LayerInterface.resize(Nz + 1);
    FormatText().Read(file_vertical_levels, LayerInterface);

    config.PeekValue("y_min", y_min);
    config.PeekValue("Delta_y", "positive", Delta_y);
    config.PeekValue("Ny", "positive", Ny);

    config.PeekValue("x_min", x_min);
    config.PeekValue("Delta_x", "positive", Delta_x);
    config.PeekValue("Nx", "positive", Nx);
  }


  //! Empty method that should check the configuration in derived classes.
  template<class T>
  void BaseModel<T>::CheckConfiguration()
  {
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Allocates memory.
  /*! Allocates the grids and the concentration Data.
   */
  template<class T>
  void BaseModel<T>::Allocate()
  {
    GridS3D = RegularGrid<T>(Ns);
    GridS4D = RegularGrid<T>(Ns);
    
    GridZ3D = RegularGrid<T>(Nz);
    GridZ3D_interf = RegularGrid<T>(Nz + 1);
    GridZ4D = RegularGrid<T>(Nz);
    GridZ4D_interf = RegularGrid<T>(Nz + 1);
    
    GridX2D = RegularGrid<T>(x_min, Delta_x, Nx);
    GridX3D = RegularGrid<T>(x_min, Delta_x, Nx);
    GridX3D_interf = RegularGrid<T>(x_min - Delta_x / 2., Delta_x, Nx + 1);
    GridX4D = RegularGrid<T>(x_min, Delta_x, Nx);
    GridX4D_interf = RegularGrid<T>(x_min - Delta_x / 2., Delta_x, Nx + 1);

    GridY2D = RegularGrid<T>(y_min, Delta_y, Ny);
    GridY3D = RegularGrid<T>(y_min, Delta_y, Ny);
    GridY3D_interf = RegularGrid<T>(y_min - Delta_y / 2., Delta_y, Ny + 1);
    GridY4D = RegularGrid<T>(y_min, Delta_y, Ny);
    GridY4D_interf = RegularGrid<T>(y_min - Delta_y / 2., Delta_y, Ny + 1);
    
    // None of the grids should be duplicated in order to save memory.
    GridS4D.SetVariable(0);
    GridS4D.SetDuplicate(false);

    GridZ3D.SetVariable(0);
    GridZ3D.SetDuplicate(false);
    GridZ3D_interf.SetVariable(0);
    GridZ3D_interf.SetDuplicate(false);
    GridZ4D.SetVariable(1);
    GridZ4D.SetDuplicate(false);
    GridZ4D_interf.SetVariable(1);
    GridZ4D_interf.SetDuplicate(false);

    GridX2D.SetVariable(1);
    GridX2D.SetDuplicate(false);
    GridX3D.SetVariable(2);
    GridX3D.SetDuplicate(false);
    GridX3D_interf.SetVariable(2);
    GridX3D_interf.SetDuplicate(false);
    GridX4D.SetVariable(3);
    GridX4D.SetDuplicate(false);
    GridX4D_interf.SetVariable(3);
    GridX4D_interf.SetDuplicate(false);

    GridY2D.SetVariable(0);
    GridY2D.SetDuplicate(false);
    GridY3D.SetVariable(1);
    GridY3D.SetDuplicate(false);
    GridY3D_interf.SetVariable(1);
    GridY3D_interf.SetDuplicate(false);
    GridY4D.SetVariable(2);
    GridY4D.SetDuplicate(false);
    GridY4D_interf.SetVariable(2);
    GridY4D_interf.SetDuplicate(false);

    Concentration.Resize(GridS4D, GridZ4D, GridY4D, GridX4D);
  }
  

  //! Model initialization.
  /*! It reads the configuration and allocates memory.
   */
  template<class T>
  void BaseModel<T>::Init()
  {
    this->ReadConfiguration();
    this->CheckConfiguration();
 
    this->Allocate();

    // Vertical layers heights.
    FormatText().Read(file_vertical_levels, GridZ3D_interf);
    for (int k=0; k<Nz; k++)
      GridZ4D(k) = GridZ3D(k)
	= (GridZ3D_interf(k) + GridZ3D_interf(k+1)) / 2.0;

    for (int k=0; k<Nz+1; k++)
      GridZ4D_interf(k) = GridZ3D_interf(k);

    BaseModel<T>::SetDate(Date_min);
  }


  /*! \brief Empty method that should initialize derived models at the
    beginning of each time step. */
  template<class T>
  void BaseModel<T>::InitStep()
  {
  }
  

  //! Moves the model to a given date.
  /*! This method prepares the model for a time integration at a given
    date. It should be called before InitStep and Forward.
    \param date date.
  */
  template<class T>
  void BaseModel<T>::SetDate(Date date)
  {
    // Sets the dates involved during the integration.
    previous_date = date;
    current_date = date;
    next_date = date;
    next_date.AddSeconds(Delta_t);
  }
  

  //! Moves back to the beginning of the previous step.
  /*! Moves back to the beginning of the previous step without changes in the
    input data and in computed concentrations.
  */
  template<class T>
  void BaseModel<T>::StepBack()
  {
    SubtractTime(Delta_t);
    this->step--;
  }


  //! Moves back to the beginning of the previous step.
  /*! Moves back to the beginning of the previous step without changes in the
    input data. Concentrations are set to \a concentration.
    \param concentration the concentration array to be set.
  */
  template<class T>
  void BaseModel<T>::StepBack(const Array<T, 4>& concentration)
  {
    this->Concentration.GetArray() = concentration;

    SubtractTime(Delta_t);
    this->step--;
  }
 

  /////////////////
  // INTEGRATION //
  /////////////////


  /*! \brief Empty methods that should perform one step forward (time
    integration over 'Delta_t') in derived models. */
  template<class T>
  void BaseModel<T>::Forward()
  {
  }

  
  //! Prepares for backward integration.
  /*! It sets flag for backward integration, then allocates memories for
    adjoint concentration data, finally it initilizes the adjoint
    concentration data to zero.
    \param flag true for backward integration; false for forward integration.
  */
  template<class T>
  void BaseModel<T>::SetBackward(bool flag)
  {
    backward = flag;
    if (backward)
      {
	Concentration_ccl.Resize(GridS4D, GridZ4D, GridY4D, GridX4D);
	Concentration_ccl.GetArray() = T(0.);
      }
  }


  /*! \brief Empty methods that should perform one step backward (time
    integration over 'Delta_t') in derived adjoint models. */
  template<class T>
  void BaseModel<T>::Backward()
  {
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////

  
  //! Returns the name of the configuration file.
  /*!
    \return The name of the configuration file.
  */
  template<class T>
  string BaseModel<T>::GetConfigurationFile() const
  {
    return file_config;
  }


  //! Is backward integration?
  /*!
    \return True for backward integration; false for forward integration.
  */
  template<class T>
  bool BaseModel<T>::IsBackward()
  {
    return backward;
  }


  //! Returns the time step.
  /*!
    \return Time step in seconds.
  */
  template<class T>
  T BaseModel<T>::GetDelta_t() const
  {
    return Delta_t;
  }


  //! Returns the number of time steps.
  /*!
    \return The number of time steps.
  */
  template<class T>
  int BaseModel<T>::GetNt() const
  {
    return Nt;
  }


  //! Returns the name of the file that describes the species.
  /*!
    \return Name of the file that describes the species.
  */
  template<class T>
  string BaseModel<T>::GetSpeciesFile() const
  {
    return file_species;
  }


  //! Returns the number of species.
  /*!
    \return The number of species.
  */
  template<class T>
  int BaseModel<T>::GetNs() const
  {
    return Ns;
  }


  //! Returns the number of species.
  /*!
    \return The number of species.
  */
  template<class T>
  int BaseModel<T>::GetNumSpecies() const
  {
    return Ns;
  }


  //! Returns the number of aerosol species.
  /*!
    \return The number of aerosol species.
  */
  template<class T>
  int BaseModel<T>::GetNs_aer() const
  {
    return Ns_aer;
  }
  

  //! Returns the number of aerosol species.
  /*!
    \return The number of aerosol species.
  */
  template<class T>
  int BaseModel<T>::GetNumSpecies_aer() const
  {
    return Ns_aer;
  }
  

  //! Returns the number of bins (aerosols).
  /*!
    \return The number of bins.
  */
  template<class T>
  int BaseModel<T>::GetNbin_aer() const
  {
    return Nbin_aer;
  }
  

  //! Returns simulation starting date.
  /*!
    \return The simulation starting date.
  */
  template<class T>
  Date BaseModel<T>::GetDate_min() const
  {
    return Date_min;
  }


  //! Returns the number of vertical layers.
  /*!
    \return The number of vertical layers.
  */
  template<class T>
  int BaseModel<T>::GetNz() const
  {
    return Nz;
  }
  

  //! Returns the domain origin along y.
  /*!
    \return The domain origin along y.
    \note In case the domain is split into cells, the origin is the center of
    the first cell.
  */
  template<class T>
  T BaseModel<T>::GetY_min() const
  {
    return y_min;
  }


  //! Returns the space step along y.
  /*!
    \return Space step along y.
  */
  template<class T>
  T BaseModel<T>::GetDelta_y() const
  {
    return Delta_y;
  }


  //! Returns the number of points along y.
  /*!
    \return The number of points along y.
  */
  template<class T>
  int BaseModel<T>::GetNy() const
  {
    return Ny;
  }


  //! Returns the domain origin along x.
  /*!
    \return The domain origin along x.
    \note In case the domain is split into cells, the origin is the center of
    the first cell.
  */
  template<class T>
  T BaseModel<T>::GetX_min() const
  {
    return x_min;
  }


  //! Returns the space step along x.
  /*!
    \return Space step along x.
  */
  template<class T>
  T BaseModel<T>::GetDelta_x() const
  {
    return Delta_x;
  }
  

  //! Returns the number of points along x.
  /*!
    \return The number of points along x.
  */
  template<class T>
  int BaseModel<T>::GetNx() const
  {
    return Nx;
  }


  //! Returns the species list.
  /*! The species order in the model (e.g. in the concentrations Data) is the
    same as the order of the returned list.
    \return The species list.
  */
  template<class T>
  vector<string> BaseModel<T>::GetSpeciesList() const
  {
    return species_list;
  }


  //! Returns the name of a species on the basis of its index.
  /*!
    \param i species index (in the model).
    \return Species name.
  */
  template<class T>
  string BaseModel<T>::GetSpeciesName(int i) const
  {
    return species_list[i];
  }


  //! Checks whether a species is a model species.
  /*!
    \param name species name.
    \return True if the species is a model species, false otherwise.
  */
  template<class T>
  bool BaseModel<T>::IsSpecies(string name) const
  {
    return find(species_list.begin(), species_list.end(), name)
      != species_list.end();
  }


  //! Returns the index of a given species.
  /*! If the species is not found, an exception is thrown.
    \param species the species name.
    \return The index (in the model) of the species named \a species.
  */
  template<class T>
  int BaseModel<T>::GetSpeciesIndex(string species) const
  {
    int index = 0;
    int length = species_list.size();
    while (index < length && species_list[index] != species)
      index++;
    if (index == length)
      throw string("Species \"") + species + "\" unknown.";
    return index;
  }
  
  
  //! Returns the index in a given vector of a species.
  /*! If the species is not found, an exception is thrown.
    \param species the species name.
    \param ref_species_list a vector of species names.
    \return The index in 'ref_species_list' of the species named \a species.
  */
  template<class T>
  int BaseModel<T>::GetSpeciesIndex(string species,
				    const vector<string>& ref_species_list)
    const
  {
    int index = 0;
    int length = ref_species_list.size();
    while (index < length && ref_species_list[index] != species)
      index++;
    if (index == length)
      throw string("Species \"") + species + "\" unknown.";
    return index;
  }


  //! Returns the index of a species in a given field.
  /*! If the species is not found, an exception is thrown.
    \param field the field name.
    \param species the species name.
    \return The index in field \a field of the species named \a species.
  */
  template<class T>
  int BaseModel<T>::GetSpeciesIndex(string field, string species)
  {
    return GetSpeciesIndex(species, GetSpeciesList(field));
  }


  //! Returns the species list of aerosol species.
  /*! The species order in the model (e.g. in the concentrations Data) is the
    same as the order of the returned list.
    \return The species list.
  */
  template<class T>
  vector<string> BaseModel<T>::GetSpeciesList_aer() const
  {
    return species_list_aer;
  }


  //! Returns the name of an aerosol on the basis of its index.
  /*!
    \param i species index (in the model).
    \return Aerosol name.
  */
  template<class T>
  string BaseModel<T>::GetSpeciesName_aer(int i) const
  {
    return species_list_aer[i];
  }


  //! Checks whether a species is a model aerosol species.
  /*!
    \param name species name.
    \return True if the species is a model aerosol species, false otherwise.
  */
  template<class T>
  bool BaseModel<T>::IsSpecies_aer(string name) const
  {
    return find(species_list_aer.begin(), species_list_aer.end(), name)
      != species_list_aer.end();
  }


  //! Returns the index of an aerosol species.
  /*! If the species is not found, an exception is thrown.
    \param species the species name.
    \return The index (in the model) of the species named \a species.
  */
  template<class T>
  int BaseModel<T>::GetSpeciesIndex_aer(string species) const
  {
    return GetSpeciesIndex(species, species_list_aer);
  }


  //! Returns the species list associated with a given field.
  /*! The species order in the model is the same as the order of the returned
    list.
    \param field the field name.
    \return The species list associated with field \a field.
  */
  template<class T>
  const vector<string>& BaseModel<T>::GetSpeciesList(string field)
  {
    map<string, vector<string>* >::const_iterator
      iter = field_species.find(field);
    if (iter == field_species.end())
      throw string("Field \"") + field + "\" is not indiced by species.";
    return *iter->second;
  }


  //! Returns the bins list associated with a given field.
  /*! The bins order in the model is the same as the order of the returned
    list.
    \param field the field name.
    \return The bins list associated with field \a field.
  */
  template<class T>
  const vector<int>& BaseModel<T>::GetBinsList_aer(string field)
  {
    map<string, vector<int>* >::const_iterator
      iter = field_bins.find(field);
    if (iter == field_bins.end())
      throw string("Field \"") + field + "\" is not indiced by bins.";
    return *iter->second;
  }


  //! Returns coordinates along x.
  /*!
    \return Coordinates along x in a 1D-array.
  */
  template<class T>
  Array<T, 1>& BaseModel<T>::GetGridXArray1D()
  {
    return GridX2D.GetArray();
  }


  //! Returns coordinates along y.
  /*!
    \return Coordinates along y in a 1D-array.
  */
  template<class T>
  Array<T, 1>& BaseModel<T>::GetGridYArray1D()
  {
    return GridY2D.GetArray();
  }


  //! Returns coordinates along z.
  /*!
    \return Coordinates along z in a 1D-array.
  */
  template<class T>
  Array<T, 1>& BaseModel<T>::GetGridZArray1D()
  {
    return GridZ3D.GetArray();
  }


  //! Returns interface altitudes.
  /*!
    \return Altitudes of layer interfaces.
  */
  template<class T>
  Array<T, 1>& BaseModel<T>::GetLayerInterface()
  {
    return LayerInterface;
  }


  //! Returns the date of the previous time-step.
  /*!
    \return The date of the previous time-step.
  */
  template<class T>
  Date BaseModel<T>::GetPreviousDate() const
  {
    return previous_date;
  }


  //! Returns the date of the current time-step.
  /*!
    \return The date of the current time-step.
  */
  template<class T>
  Date BaseModel<T>::GetCurrentDate() const
  {
    return current_date;
  }


  /*! \brief Returns the number of seconds between the simulation beginning
    and the current date. */
  /*!
    \return The number of seconds between the simulation beginning and the
    current date.
  */
  template<class T>
  T BaseModel<T>::GetCurrentTime() const
  {
    return current_time;
  }


  //! Returns the date of the next time-step.
  /*!
    \return The date of the next time-step.
  */
  template<class T>
  Date BaseModel<T>::GetNextDate() const
  {
    return next_date;
  }


  //! Checks whether a field is part of the model interface.
  /*!
    \param field field name.
    \return True if the field is part of the interface, false otherwise.
  */
  template<class T>
  bool BaseModel<T>::HasField(string field)
  {
    return D2_map.find(field) != D2_map.end()
      || A2_map.find(field) != A2_map.end()
      || D3_map.find(field) != D3_map.end()
      || A3_map.find(field) != A3_map.end()
      || D4_map.find(field) != D4_map.end()
      || A4_map.find(field) != A4_map.end()
      || D5_map.find(field) != D5_map.end()
      || A5_map.find(field) != A5_map.end();
  }


  /*! \brief Returns information about the availability of a field in the
    model interface. */
  /*! If the field is found in the model interface, a string is returned with
    the description of the field in the model: dimension and availability as
    an array or as a Data instance. Otherwise, the returned string says that
    the field is not part of the interface.
    \param field field name.
    \return A string that says whether the field is available, and in which
    form (array or Data, along with its dimension) if it is part of the
    interface.
  */
  template<class T>
  string BaseModel<T>::FindField(string field)
  {
    vector<string> type;

    // Is it a 2D field?
    if (D2_map.find(field) != D2_map.end())
      type.push_back("2D data");
    else if (A2_map.find(field) != A2_map.end())
      type.push_back("2D array");
    
    // Is it a 3D field?
    if (D3_map.find(field) != D3_map.end())
      type.push_back("3D data");
    else if (A3_map.find(field) != A3_map.end())
      type.push_back("3D array");

    // Is it a 4D field?
    if (D4_map.find(field) != D4_map.end())
      type.push_back("4D data");
    else if (A4_map.find(field) != A4_map.end())
      type.push_back("4D array");

    // Is it a 5D field?
    if (D5_map.find(field) != D5_map.end())
      type.push_back("5D data");
    else if (A5_map.find(field) != A5_map.end())
      type.push_back("5D array");

    // The field was not found.
    if (type.size() == 0)
      return string("Field \"") + field
	+ string("\" is not part of the interface.");
    
    // The field was found.
    string message = string("Field \"") + field
      + string("\" may be found in ") + type[0];
    if (type.size() == 1)
      return message + ".";
    for (unsigned int i = 1; i < type.size() - 1; i++)
      message += string(", ") + type[i];
    return message + string(" and ") + type[type.size() - 1] + ".";
  }


  //! Returns the number of dimensions of a field in the model interface.
  /*! If the field is found in the model interface, its dimension is returned.
    Otherwise it raises an exception.
    \param field field name.
    \return An integer that indicates the field dimension (array or Data): 2
    for 2D arrays, 3 for 3D arrays, ...
  */
  template<class T>
  int BaseModel<T>::GetFieldDimension(string field)
  {
    // Is it a 2D field?
    if (D2_map.find(field) != D2_map.end())
      return 2;
    else if (A2_map.find(field) != A2_map.end())
      return 2;
    
    // Is it a 3D field?
    if (D3_map.find(field) != D3_map.end())
      return 3;
    else if (A3_map.find(field) != A3_map.end())
      return 3;

    // Is it a 4D field?
    if (D4_map.find(field) != D4_map.end())
      return 4;
    else if (A4_map.find(field) != A4_map.end())
      return 4;

    // Is it a 5D field?
    if (D5_map.find(field) != D5_map.end())
      return 5;
    else if (A5_map.find(field) != A5_map.end())
      return 5;

    throw string("Unable to find field \"") + field
      + "\" in the model interface.";
  }


  //! Returns the array that stores a 2D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 2D array that stores \a field.
  */
  template<class T>
  Array<T, 2>& BaseModel<T>::A2(string field)
  {
    if (A2_map.find(field) == A2_map.end())
      if (D2_map.find(field) == D2_map.end())
	throw string("Field \"") + field + string("\" cannot be found in 2D")
	  + string(" arrays. ") + FindField(field);
      else
	return D2_map[field]->GetArray();
    else
      return *A2_map[field];
  }

  
  //! Returns the Data instance that stores a 2D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 2D Data instance that stores \a field.
  */
  template<class T>
  Data<T, 2>& BaseModel<T>::D2(string field)
  {
    if (D2_map.find(field) == D2_map.end())
      throw string("Field \"") + field + string("\" cannot be found in 2D")
	+ string(" data. ") + FindField(field);
    return *D2_map[field];
  }

  
  //! Returns the array that stores a 3D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 3D array that stores \a field.
  */
  template<class T>
  Array<T, 3>& BaseModel<T>::A3(string field)
  {
    if (A3_map.find(field) == A3_map.end())
      if (D3_map.find(field) == D3_map.end())
	throw string("Field \"") + field + string("\" cannot be found in 3D")
	  + string(" arrays. ") + FindField(field);
      else
	return D3_map[field]->GetArray();
    else
      return *A3_map[field];
  }

  
  //! Returns the Data instance that stores a 3D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 3D Data instance that stores \a field.
  */
  template<class T>
  Data<T, 3>& BaseModel<T>::D3(string field)
  {
    if (D3_map.find(field) == D3_map.end())
      throw string("Field \"") + field + string("\" cannot be found in 3D")
	+ string(" data. ") + FindField(field);
    return *D3_map[field];
  }

  
  //! Returns the array that stores a 4D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 4D array that stores \a field.
  */
  template<class T>
  Array<T, 4>& BaseModel<T>::A4(string field)
  {
    if (A4_map.find(field) == A4_map.end())
      if (D4_map.find(field) == D4_map.end())
	throw string("Field \"") + field + string("\" cannot be found in 4D")
	  + string(" arrays. ") + FindField(field);
      else
	return D4_map[field]->GetArray();
    else
      return *A4_map[field];
  }

  
  //! Returns the Data instance that stores a 4D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 4D Data instance that stores \a field.
  */
  template<class T>
  Data<T, 4>& BaseModel<T>::D4(string field)
  {
    if (D4_map.find(field) == D4_map.end())
      throw string("Field \"") + field + string("\" cannot be found in 4D")
	+ string(" data. ") + FindField(field);
    return *D4_map[field];
  }


  //! Returns the array that stores a 5D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 5D array that stores \a field.
  */
  template<class T>
  Array<T, 5>& BaseModel<T>::A5(string field)
  {
    if (A5_map.find(field) == A5_map.end())
      if (D5_map.find(field) == D5_map.end())
	throw string("Field \"") + field + string("\" cannot be found in 5D")
	  + string(" arrays. ") + FindField(field);
      else
	return D5_map[field]->GetArray();
    else
      return *A5_map[field];
  }

  
  //! Returns the Data instance that stores a 5D field.
  /*! If the field is not found, an exception is thrown.
    \param field field name.
    \return The 5D Data instance that stores \a field.
  */
  template<class T>
  Data<T, 5>& BaseModel<T>::D5(string field)
  {
    if (D5_map.find(field) == D5_map.end())
      throw string("Field \"") + field + string("\" cannot be found in 5D")
	+ string(" data. ") + FindField(field);
    return *D5_map[field];
  }


  //! Returns the concentrations Data.
  /*!
    \return The concentrations Data.
  */
  template<class T>
  Data<T, 4>& BaseModel<T>::GetConcentration()
  {
    return Concentration;
  }


  //! Returns the concentrations Data.
  /*!
    \return The concentrations Data.
  */
  template<class T>
  const Data<T, 4>& BaseModel<T>::GetConcentration() const
  {
    return Concentration;
  }


  //! Returns the adjoint concentrations Data.
  /*!
    \return The adjoint concentrations Data.
  */
  template<class T>
  Data<T, 4>& BaseModel<T>::GetConcentration_ccl()
  {
    return Concentration_ccl;
  }


  //! Returns the adjoint concentrations Data.
  /*!
    \return The adjoint concentrations Data.
  */
  template<class T>
  const Data<T, 4>& BaseModel<T>::GetConcentration_ccl() const
  {
    return Concentration_ccl;
  }


  //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  Data<T, 5>& BaseModel<T>::GetConcentration_aer()
  {
    return Concentration_aer;
  }


  //! Returns the concentrations Data for aerosols.
  /*!
    \return The concentrations Data for aerosols.
  */
  template<class T>
  const Data<T, 5>& BaseModel<T>::GetConcentration_aer() const
  {
    return Concentration_aer;
  }


  //! Computes mean concentration over a given volume (virtual method).
  /*!
    \return The concentration averaged over a given volume.
  */
  template<class T>
  T BaseModel<T>::GetIntegratedConcentration(int species, T z, T y, T x,
					     T lz, T ly, T lx)
  {
    throw string("\"BaseModel<T>::GetIntegratedConcentration(int, T, T, T, ")
      + "T, T, T)\" is not defined.";
  }


  //! Computes concentration at a given point (virtual method).
  /*!
    \return The concentration at the point.
  */
  template<class T>
  T BaseModel<T>::GetConcentration(int species, T z, T y, T x)
  {
    throw string("\"BaseModel<T>::GetConcentration(int, T, T, T)\"")
      + " is not defined.";
  }


  //! Computes concentration at a given point (virtual method) for aerosols.
  /*!
    \return The aerosol concentration at the point.
  */
  template<class T>
  T BaseModel<T>::GetConcentration_aer(int species, int diameter,
				       T z, T y, T x)
  {
    throw string("\"BaseModel<T>::GetConcentration_aer(int, int, T, T, T)\"")
      + " is not defined.";
  }


  //! Computes concentration for a list of species and levels (virtual
  //! method).
  template<class T>
  void BaseModel<T>::ComputeConcentration(const vector<int>& species_index,
					  const vector<int>& levels)
  {
  }


  //! Computes concentration on the whole grid of the model (virtual method).
  template<class T>
  void BaseModel<T>::ComputeConcentration()
  {
  }


  //! Checks whether a field is managed by the model.
  /*! \note A field managed by the model is a field that is read or computed
    by the model.
    \param field field name.
    \return True if the field is managed by the model, false otherwise.
  */
  template<class T>
  bool BaseModel<T>::IsManaged(string field) const
  {
    map<string, bool>::const_iterator iter = option_manage.find(field);
    if (iter == option_manage.end())
      throw string("Field \"") + field + "\" is unknown.";
    return iter->second;
  }


  //! Specifies whether a given field should be managed by the model.
  /*! \note A field managed by the model is a field that is read or computed
    by the model.
    \param field field name.
    \param status true if the field should be managed by the model, false
    otherwise.
  */
  template<class T>
  void BaseModel<T>::Manage(string field, bool status)
  {
    map<string, bool>::iterator iter = option_manage.find(field);
    if (iter == option_manage.end())
      throw string("Field \"") + field + "\" is unknown.";
    iter->second = status;
  }


  //! Sets date of the current time-step to a given date.
  /*!
    \param date date.
    \note Fields are not synchronized.
  */
  template<class T>
  void BaseModel<T>::SetCurrentDate(Date date)
  {
    // Sets the dates involved during the integration.
    BaseModel<T>::SetDate(date);
  }


  //! Adds time to current time.
  /*! Sets the previous date to the (old) current date, adds \a time seconds
    to the current date, and sets the next date to the (new) current date plus
    \a time seconds.
    \param time seconds to be added.
  */
  template<class T>
  void BaseModel<T>::AddTime(T time)
  {
    previous_date = current_date;

    current_date.AddSeconds(time);

    next_date = current_date;
    next_date.AddSeconds(time);

    current_time += time;
  }
  

  //! Subtracts time to current time.
  /*! Sets the next date to the (old) current date, subtracts \a time seconds
    to the current date, and sets the previous date to the (new) current date
    minus \a time seconds.
    \param time seconds to be subtracted.
  */
  template<class T>
  void BaseModel<T>::SubtractTime(T time)
  {
    next_date = current_date;

    current_date.AddSeconds(-time);

    previous_date = current_date;
    previous_date.AddSeconds(-time);

    current_time -= time;
  }
  

  ///////////////////////
  // PROTECTED METHODS //
  ///////////////////////


  //! Initializes data at a given date.
  /*! Reads, on file, data associated with \a field in group \a section. It
    reads the two steps surrounding the date \a date. Field \a CurrentData is
    then computed by linear interpolation or set to \a FileData_f, depending
    on \a interpolated.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (output) data stored on file at the step before \a date.
    \param FileData_f (output) data stored on file at the step after \a date.
    \param date date at which data are computed.
    \param CurrentData (output) data at date \a date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \note The input file is 'this->input_files[section](field)'.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::InitData(string section, string field,
			      Data<T, N>& FileData_i, Data<T, N>& FileData_f,
			      Date date, Data<T, N>& CurrentData,
			      bool interpolated)
  {
    InitData(input_files[section](field), input_files[section].GetDateMin(),
	     input_files[section].GetDelta_t(), FileData_i, FileData_f,
	     date, CurrentData, interpolated);
  }
  
  
  //! Initializes sub-data at a given date.
  /*! Reads, on file, data associated with \a field in group \a section. It
    reads the two steps surrounding the date \a date. Field \a CurrentData is
    then computed by linear interpolation or set to \a FileData_f, depending
    on \a interpolated. Output data is stored in the ith sub-array of output
    Data instances, where i is the index along the first dimension.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (output) data stored on file at the step before \a date.
    \param FileData_f (output) data stored on file at the step after \a date.
    \param date date at which data are computed.
    \param i index, along the first dimension, where output data is stored.
    \param CurrentData (output) data at date \a date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \note The input file is 'this->input_files[section](field)'.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::InitData(string section, string field,
			      Data<T, N>& FileData_i, Data<T, N>& FileData_f,
			      Date date, int index, Data<T, N>& CurrentData,
			      bool interpolated)
  {
    InitData(input_files[section](field), input_files[section].GetDateMin(),
	     input_files[section].GetDelta_t(), FileData_i, FileData_f, date,
	     index, CurrentData, interpolated);
  }
  
  
  //! Initializes sub-data at a given date.
  /*! Reads, on file, data associated with \a field in group \a section. It
    reads the two steps surrounding the date \a date and it interpolates
    linearly to compute the data at date \a date. Output data is stored in the
    sub-array at (first_index, second_index) of output Data instances, where
    \a first_index is the index along the first dimension, and \a second_index
    is the index along the second dimension.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (output) data stored on file at the step before \a date.
    \param FileData_f (output) data stored on file at the step after \a date.
    \param date date at which data are computed.
    \param first_index index, along the first dimension, where output data is
    stored.
    \param second_index index, along the second dimension, where output data
    is stored.
    \param CurrentData (output) data at date \a date.
    \note The input file is 'this->input_files[section](field)'.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::InitData(string section, string field,
			      Data<T, N>& FileData_i, Data<T, N>& FileData_f,
			      Date date, int first_index, int second_index,
			      Data<T, N>& CurrentData)
  {
    InitData(input_files[section](field), input_files[section].GetDateMin(),
	     input_files[section].GetDelta_t(), FileData_i, FileData_f, date,
	     first_index, second_index, CurrentData);
  }
  
  
  //! Initializes data at a given date.
  /*! Reads data in file \a input_file. It reads the two steps surrounding the
    date \a date. Field \a CurrentData is then computed by linear
    interpolation or set to \a FileData_f, depending on \a interpolated.
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (output) data stored on file at the step before \a date.
    \param FileData_f (output) data stored on file at the step after \a date.
    \param date date at which data are computed.
    \param CurrentData (output) data at date \a date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \warning If the input-file name is a number, the data is set to this
    number.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::InitData(string input_file, Date date_min_file,
			      T Delta_t_file, Data<T, N>& FileData_i,
			      Data<T, N>& FileData_f, Date date,
			      Data<T, N>& CurrentData,
			      bool interpolated)
  {
    if (is_num(input_file))
      {
	T value = to_num<T>(input_file);
	FileData_i.Fill(value);
	FileData_f.Fill(value);
	CurrentData.Fill(value);
	return;
      }

    T time_distance = date.GetSecondsFrom(date_min_file);
    int record = int(floor(time_distance / Delta_t_file));
    FormatBinary<float>().ReadRecord(input_file, record, FileData_i);
    FormatBinary<float>().ReadRecord(input_file, record + 1, FileData_f);

    if (interpolated)
      Interpolate(date_min_file, Delta_t_file, record, FileData_i, FileData_f,
		  date, CurrentData);
    else
      CurrentData.GetArray() = FileData_f.GetArray();
  }
  
  
  //! Initializes sub-data at a given date.
  /*! Reads data in file \a input_file. It reads the two steps surrounding the
    date \a date. Field \a CurrentData is then computed by linear
    interpolation or set to \a FileData_f, depending on \a interpolated.
    Output data is stored in the ith sub-array of output Data instances, where
    i is the index along the first dimension.
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (output) data stored on file at the step before \a date.
    \param FileData_f (output) data stored on file at the step after \a date.
    \param date date at which data are computed.
    \param i index, along the first dimension, where output data is stored.
    \param CurrentData (output) data at date \a date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \warning If the input-file name is a number, the data is set to this
    number.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::InitData(string input_file, Date date_min_file,
			      T Delta_t_file, Data<T, N>& FileData_i,
			      Data<T, N>& FileData_f, Date date, int index,
			      Data<T, N>& CurrentData,
			      bool interpolated)
  {
    // Extracts sub-data first.
    TinyVector<int, N-1> new_shape;
    for (int i = 0; i < N - 1; i++)
      new_shape(i) = FileData_i.GetArray().shape()(i + 1);
    unsigned int position = index * FileData_i.GetNbElements()
      / FileData_i.GetLength(0);
    Data<T, N-1> FileData_extract_i(&FileData_i.GetData()[position],
				    new_shape);
    Data<T, N-1> FileData_extract_f(&FileData_f.GetData()[position],
				    new_shape);
    Data<T, N-1> CurrentData_extract(&CurrentData.GetData()[position],
				     new_shape);

    InitData(input_file, date_min_file, Delta_t_file, FileData_extract_i,
	     FileData_extract_f, date, CurrentData_extract, interpolated);
  }


  //! Initializes sub-data at a given date.
  /*! Reads data in file \a input_file. It reads the two steps surrounding the
    date \a date and it interpolates linearly to compute the data at date \a
    date. Output data is stored in the sub-array at (first_index,
    second_index) of output Data instances, where \a first_index is the index
    along the first dimension, and \a second_index is the index along the
    second dimension.
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (output) data stored on file at the step before \a date.
    \param FileData_f (output) data stored on file at the step after \a date.
    \param date date at which data are computed.
    \param first_index index, along the first dimension, where output data is
    stored.
    \param second_index index, along the second dimension, where output data
    is stored.
    \param CurrentData (output) data at date \a date.
    \warning If the input-file name is a number, the data is set to this
    number.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::InitData(string input_file, Date date_min_file,
			      T Delta_t_file, Data<T, N>& FileData_i,
			      Data<T, N>& FileData_f, Date date,
			      int first_index, int second_index,
			      Data<T, N>& CurrentData)
  {
    // Extracts sub-data first.
    TinyVector<int, N-2> new_shape;
    for (int i = 0; i < N - 2; i++)
      new_shape(i) = FileData_i.GetArray().shape()(i + 2);
    unsigned int position = first_index * FileData_i.GetNbElements()
      + second_index * FileData_i.GetNbElements() / FileData_i.GetLength(1);
    position /= FileData_i.GetLength(0);
    Data<T, N-2> FileData_extract_i(&FileData_i.GetData()[position],
				    new_shape);
    Data<T, N-2> FileData_extract_f(&FileData_f.GetData()[position],
				    new_shape);
    Data<T, N-2> CurrentData_extract(&CurrentData.GetData()[position],
				     new_shape);

    InitData(input_file, date_min_file, Delta_t_file, FileData_extract_i,
	     FileData_extract_f, date, CurrentData_extract);
  }


  //! Updates data at current date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that surround the previous date. If they do not surround the current
    date, they are updated (i.e. read on file). The input file is associated
    with \a field in group \a section.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param CurrentData (output) data at current date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string section, string field,
				Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				Data<T, N>& CurrentData,
				bool interpolated)
  {
    UpdateData(input_files[section](field), input_files[section].GetDateMin(),
	       input_files[section].GetDelta_t(), FileData_i, FileData_f,
	       CurrentData, interpolated);
  }


  //! Updates sub-data at current date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that surround the previous date. If they do not surround the current
    date, they are updated (i.e. read on file). The input file is associated
    with \a field in group \a section. Output data is stored in the ith
    sub-array of output Data instances, where i is the index along the first
    dimension.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param i index, along the first dimension, where output data is stored.
    \param CurrentData (output) data at current date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string section, string field,
				Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				int index, Data<T, N>& CurrentData,
				bool interpolated)
  {
    UpdateData(input_files[section](field), input_files[section].GetDateMin(),
	       input_files[section].GetDelta_t(), FileData_i, FileData_f,
	       index, CurrentData, interpolated);
  }


  //! Updates sub-data at current date.
  /*! Interpolates data at the current date, on the basis of \a FileData_i and
    \a FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read
    on file) that surround the previous date. If they do not surround the
    current date, they are updated (i.e. read on file). The input file is
    associated with \a field in group \a section. Output data is stored in the
    sub-array at (first_index, second_index) of output Data instances, where
    \a first_index is the index along the first dimension, and \a second_index
    is the index along the second dimension.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param first_index index, along the first dimension, where output data is
    stored.
    \param second_index index, along the second dimension, where output data
    is stored.
    \param CurrentData (output) data at current date.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string section, string field,
				Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				int first_index, int second_index,
				Data<T, N>& CurrentData)
  {
    UpdateData(input_files[section](field), input_files[section].GetDateMin(),
	       input_files[section].GetDelta_t(), FileData_i, FileData_f,
	       first_index, second_index, CurrentData);
  }


  //! Updates data at current date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that surround the previous date. If they do not surround the current
    date, they are updated (i.e. read on file).
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param CurrentData (output) data at current date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date.
    \warning If the input-file name is a number, the data is set to this
    number.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string input_file, Date date_min_file,
				T Delta_t_file, Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				Data<T, N>& CurrentData,
				bool interpolated)
  {
    if (is_num(input_file))
      {
	T value = to_num<T>(input_file);
	FileData_i.Fill(value);
	FileData_f.Fill(value);
	CurrentData.Fill(value);
	return;
      }

    // Computes the record of the previous step, to check whether
    // any update of FileData_{i,f} is needed.
    T time_distance = previous_date.GetSecondsFrom(date_min_file);
    int previous_record = int(floor(time_distance / Delta_t_file));
    // Computes the record to be read in the file for the current date.
    time_distance = current_date.GetSecondsFrom(date_min_file);
    int record = int(floor(time_distance / Delta_t_file));

    if (previous_record == record - 1)
      // FileData_{i,f} should be updated, and the new value of FileData_i
      // is the old value of FileData_f.
      {
	FileData_i.GetArray() = FileData_f.GetArray();
	FormatBinary<float>().ReadRecord(input_file, record + 1, FileData_f);
      }
    else if (previous_record < record)
      // All should be read.
      {
	FormatBinary<float>().ReadRecord(input_file, record, FileData_i);
	FormatBinary<float>().ReadRecord(input_file, record + 1, FileData_f);
      }

    if (interpolated)
      Interpolate(date_min_file, Delta_t_file, record, FileData_i, FileData_f,
		  current_date, CurrentData);
    else
      CurrentData.GetArray() = FileData_f.GetArray();
  }
  
  
  //! Updates sub-data at current date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that surround the previous date. If they do not surround the current
    date, they are updated (i.e. read on file). Output data is stored in the
    ith sub-array of output Data instances, where i is the index along the
    first dimension.
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param i index, along the first dimension, where output data is stored.
    \param CurrentData (output) data at current date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date.
    \warning If the input-file name is a number, the data is set to this
    number.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string input_file, Date date_min_file,
				T Delta_t_file, Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				int index, Data<T, N>& CurrentData,
				bool interpolated)
  {
    // Extracts sub-data first.
    TinyVector<int, N-1> new_shape;
    for (int i = 0; i < N - 1; i++)
      new_shape(i) = FileData_i.GetArray().shape()(i + 1);
    unsigned int position = index * FileData_i.GetNbElements()
      / FileData_i.GetLength(0);
    Data<T, N-1> FileData_extract_i(&FileData_i.GetData()[position],
				    new_shape);
    Data<T, N-1> FileData_extract_f(&FileData_f.GetData()[position],
				    new_shape);
    Data<T, N-1> CurrentData_extract(&CurrentData.GetData()[position],
				     new_shape);

    UpdateData(input_file, date_min_file, Delta_t_file, FileData_extract_i,
	       FileData_extract_f, CurrentData_extract, interpolated);
  }
  
  
  //! Updates sub-data at current date.
  /*! Interpolates data at the current date, on the basis of \a FileData_i and
    \a FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read
    on file) that surround the previous date. If they do not surround the
    current date, they are updated (i.e. read on file). Output data is stored
    in the sub-array at (first_index, second_index) of output Data instances,
    where \a first_index is the index along the first dimension, and \a
    second_index is the index along the second dimension.
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param first_index index, along the first dimension, where output data is
    stored.
    \param second_index index, along the second dimension, where output data
    is stored.
    \param CurrentData (output) data at current date.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date.
    \warning If the input-file name is a number, the data is set to this
    number.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string input_file, Date date_min_file,
				T Delta_t_file, Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				int first_index, int second_index,
				Data<T, N>& CurrentData)
  {
    // Extracts sub-data first.
    TinyVector<int, N-2> new_shape;
    for (int i = 0; i < N - 2; i++)
      new_shape(i) = FileData_i.GetArray().shape()(i + 2);
    unsigned int position = first_index * FileData_i.GetNbElements()
      + second_index * FileData_i.GetNbElements() / FileData_i.GetLength(1);
    position /= FileData_i.GetLength(0);
    Data<T, N-2> FileData_extract_i(&FileData_i.GetData()[position],
				    new_shape);
    Data<T, N-2> FileData_extract_f(&FileData_f.GetData()[position],
				    new_shape);
    Data<T, N-2> CurrentData_extract(&CurrentData.GetData()[position],
				     new_shape);

    UpdateData(input_file, date_min_file, Delta_t_file, FileData_extract_i,
	       FileData_extract_f, CurrentData_extract);
  }
  
  
  //! Updates data at current date and at next date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that surround the current date. If they do not surround the next
    date, they are updated (i.e. read on file). The input file is associated
    with \a field in group \a section.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param CurrentData_i (output) data at current date.
    \param CurrentData_f (input/output) on entry, data at current date; on
    exit, data at next date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date. On entry, \a CurrentData_f
    must be the data at current date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string section, string field,
				Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				Data<T, N>& CurrentData_i,
				Data<T, N>& CurrentData_f,
				bool interpolated)
  {
    UpdateData(input_files[section](field), input_files[section].GetDateMin(),
	       input_files[section].GetDelta_t(), FileData_i, FileData_f,
	       CurrentData_i, CurrentData_f, interpolated);
  }
  

  //! Updates sub-data at current date and at next date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that surround the current date. If they do not surround the next
    date, they are updated (i.e. read on file). The input file is associated
    with \a field in group \a section. Output data is stored in the ith
    sub-array of output Data instances, where i is the index along the first
    dimension.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param i index, along the first dimension, where output data is stored.
    \param CurrentData_i (output) data at current date.
    \param CurrentData_f (input/output) on entry, data at current date; on
    exit, data at next date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date. On entry, \a CurrentData_f
    must be the data at current date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string section, string field,
				Data<T, N>& FileData_i,
				Data<T, N>& FileData_f, int index,
				Data<T, N>& CurrentData_i,
				Data<T, N>& CurrentData_f,
				bool interpolated)
  {
    UpdateData(input_files[section](field), input_files[section].GetDateMin(),
	       input_files[section].GetDelta_t(), FileData_i, FileData_f,
	       index, CurrentData_i, CurrentData_f, interpolated);
  }
  

  //! Updates sub-data at current date and at next date.
  /*! Interpolates data at the current date and at the next date, on the basis
    of \a FileData_i and \a FileData_f. On entry, \a FileData_i and \a
    FileData_f are arrays (read on file) that surround the current date. If
    they do not surround the next date, they are updated (i.e. read on
    file). The input file is associated with \a field in group \a
    section. Output data is stored in the sub-array at (first_index,
    second_index) of output Data instances, where \a first_index is the index
    along the first dimension, and \a second_index is the index along the
    second dimension.
    \param section group to which the field belongs.
    \param field field name.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param first_index index, along the first dimension, where output data is
    stored.
    \param second_index index, along the second dimension, where output data
    is stored.
    \param CurrentData_i (output) data at current date.
    \param CurrentData_f (input/output) on entry, data at current date; on
    exit, data at next date.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date. On entry, \a CurrentData_f
    must be the data at current date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string section, string field,
				Data<T, N>& FileData_i,
				Data<T, N>& FileData_f, int first_index,
				int second_index, Data<T, N>& CurrentData_i,
				Data<T, N>& CurrentData_f)
  {
    UpdateData(input_files[section](field), input_files[section].GetDateMin(),
	       input_files[section].GetDelta_t(), FileData_i, FileData_f,
	       first_index, second_index, CurrentData_i, CurrentData_f);
  }
  

  //! Updates data at current date and at next date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read in
    file \a input_file) that surround the current date. If they do not
    surround the next date, they are updated (i.e. read in file \a
    input_file).
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param CurrentData_i (output) data at current date.
    \param CurrentData_f (input/output) on entry, data at current date; on
    exit, data at next date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the current date. On entry, \a CurrentData_f
    must be the data at current date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string input_file, Date date_min_file,
				T Delta_t_file, Data<T, N>& FileData_i,
				Data<T, N>& FileData_f,
				Data<T, N>& CurrentData_i,
				Data<T, N>& CurrentData_f,
				bool interpolated)
  {

    /*** CurrentData_i ***/

    // On entry, 'CurrentData_f' must be the data at current date.
    CurrentData_i.GetArray() = CurrentData_f.GetArray();

    /*** CurrentData_f ***/

    if (is_num(input_file))
      {
	CurrentData_f.Fill(to_num<T>(input_file));
	return;
      }

    // Computes the record of the previous step, to check whether
    // any update of FileData_{i,f} is needed. Currently, FileData_{i,f}
    // are assumed to be suited for interpolation at current_date.
    T time_distance = current_date.GetSecondsFrom(date_min_file);
    int previous_record = int(floor(time_distance / Delta_t_file));
    // Computes the record to be read in the file for the next date.
    time_distance = next_date.GetSecondsFrom(date_min_file);
    int record = int(floor(time_distance / Delta_t_file));

    if (previous_record == record - 1)
      // FileData_{i,f} should be updated, and the new value of FileData_i
      // is the old value of FileData_f.
      {
	FileData_i.GetArray() = FileData_f.GetArray();
	FormatBinary<float>().ReadRecord(input_file, record + 1, FileData_f);
      }
    else if (previous_record < record)
      // All should be read.
      {
	FormatBinary<float>().ReadRecord(input_file, record, FileData_i);
	FormatBinary<float>().ReadRecord(input_file, record + 1, FileData_f);
      }

    if (interpolated)
      Interpolate(date_min_file, Delta_t_file, record, FileData_i, FileData_f,
		  next_date, CurrentData_f);
    else
      CurrentData_f.GetArray() = FileData_f.GetArray();
  }
  
  
  //! Updates sub-data at current date and at next date.
  /*! Depending on \a interpolated, it interpolates data at the current date,
    on the basis of \a FileData_i and \a FileData_f or sets it to \a
    FileData_f. On entry, \a FileData_i and \a FileData_f are arrays (read in
    file \a input_file) that surround the current date. If they do not
    surround the next date, they are updated (i.e. read in file \a
    input_file). Output data is stored in the ith sub-array of output Data
    instances, where i is the index along the first dimension.
    \param input_file input file.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param i index, along the first dimension, where output data is stored.
    \param CurrentData_i (output) data at current date.
    \param CurrentData_f (input/output) on entry, data at current date; on
    exit, data at next date.
    \param interpolated (optional) true if \a CurrentData should be computed
    with time interpolation, false if \a CurrentData should be set to \a
    FileData_f. Default: true.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date. On entry, \a CurrentData_f
    must be the data at current date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string input_file, Date date_min_file,
				T Delta_t_file, Data<T, N>& FileData_i,
				Data<T, N>& FileData_f, int index,
				Data<T, N>& CurrentData_i,
				Data<T, N>& CurrentData_f,
				bool interpolated)
  {
    // Extracts sub-data first.
    TinyVector<int, N-1> new_shape;
    for (int i = 0; i < N - 1; i++)
      new_shape(i) = FileData_i.GetArray().shape()(i + 1);
    unsigned int position = index * FileData_i.GetNbElements()
      / FileData_i.GetLength(0);
    Data<T, N-1> FileData_extract_i(&FileData_i.GetData()[position],
				    new_shape);
    Data<T, N-1> FileData_extract_f(&FileData_f.GetData()[position],
				    new_shape);
    Data<T, N-1> CurrentData_extract_i(&CurrentData_i.GetData()[position],
				       new_shape);
    Data<T, N-1> CurrentData_extract_f(&CurrentData_f.GetData()[position],
				       new_shape);

    UpdateData(input_file, date_min_file, Delta_t_file,
	       FileData_extract_i, FileData_extract_f,
	       CurrentData_extract_i, CurrentData_extract_f, interpolated);
  }


  //! Updates sub-data at current date and at next date.
  /*! Interpolates data at the current date and at the next date, on the basis
    of \a FileData_i and \a FileData_f. On entry, \a FileData_i and \a
    FileData_f are arrays (read in file \a input_file) that surround the
    current date. If they do not surround the next date, they are updated
    (i.e. read in file \a input_file). Output data is stored in the sub-array
    at (first_index, second_index) of output Data instances, where \a
    first_index is the index along the first dimension, and \a second_index is
    the index along the second dimension.
    \param date_min_file date of the first step stored in the input file.
    \param Delta_t_file time step of the input file, in seconds.
    \param FileData_i (input/output) data stored on file at the step before
    \a date.
    \param FileData_f (input/output) data stored on file at the step after
    \a date.
    \param first_index index, along the first dimension, where output data is
    stored.
    \param second_index index, along the second dimension, where output data
    is stored.
    \param CurrentData_i (output) data at current date.
    \param CurrentData_f (input/output) on entry, data at current date; on
    exit, data at next date.
    \note The input file is 'this->input_files[section](field)'.
    \warning On entry, \a FileData_i and \a FileData_f are arrays (read on
    file) that must surround the previous date. On entry, \a CurrentData_f
    must be the data at current date.
  */
  template<class T>
  template<int N>
  void BaseModel<T>::UpdateData(string input_file, Date date_min_file,
				T Delta_t_file, Data<T, N>& FileData_i,
				Data<T, N>& FileData_f, int first_index,
				int second_index, Data<T, N>& CurrentData_i,
				Data<T, N>& CurrentData_f)
  {
    // Extracts sub-data first.
    TinyVector<int, N-2> new_shape;
    for (int i = 0; i < N - 2; i++)
      new_shape(i) = FileData_i.GetArray().shape()(i + 2);
    unsigned int position = first_index * FileData_i.GetNbElements()
      + second_index * FileData_i.GetNbElements() / FileData_i.GetLength(1);
    position /= FileData_i.GetLength(0);
    Data<T, N-2> FileData_extract_i(&FileData_i.GetData()[position],
				    new_shape);
    Data<T, N-2> FileData_extract_f(&FileData_f.GetData()[position],
				    new_shape);
    Data<T, N-2> CurrentData_extract_i(&CurrentData_i.GetData()[position],
				       new_shape);
    Data<T, N-2> CurrentData_extract_f(&CurrentData_f.GetData()[position],
				       new_shape);

    UpdateData(input_file, date_min_file, Delta_t_file,
	       FileData_extract_i, FileData_extract_f,
	       CurrentData_extract_i, CurrentData_extract_f);
  }


  //! Interpolates linearly in time.
  /*! Data on which the interpolation is based is assumed to come from a file.
    \param date_min_file starting date of the file.
    \param Delta_t_file file time-step.
    \param record record of \a FileData_i in the file.
    \param FileData_i data read in the file at record \a record.
    \param FileData_f data read in the file at record \a record + 1.
    \param date date at which data is interpolated.
    \param InterpolatedData (output) data linearly interpolated between
    \a FileData_i and \a FileData_f.
  */
  template<class T>
  void BaseModel<T>::Interpolate(Date date_min_file, T Delta_t_file,
				 int record, const Data<T, 3>& FileData_i,
				 const Data<T, 3>& FileData_f, Date date,
				 Data<T, 3>& InterpolatedData)
  {
    int k, j, i;
    int Nz = FileData_i.GetLength(0);
    int Ny = FileData_i.GetLength(1);
    int Nx = FileData_i.GetLength(2);
    
    T time_distance = date.GetSecondsFrom(date_min_file);
    T weight_f = (time_distance - T(record) * Delta_t_file) / Delta_t_file;
    T weight_i = 1. - weight_f;

    for (k = 0; k < Nz; k++)
      for (j = 0; j < Ny; j++)
	for (i = 0; i < Nx; i++)
	  InterpolatedData(k, j, i) = FileData_i(k, j, i) * weight_i
	    + FileData_f(k, j, i) * weight_f;
  }


  //! Interpolates linearly in time.
  /*! Data on which the interpolation is based is assumed to come from a file.
    \param date_min_file starting date of the file.
    \param Delta_t_file file time-step.
    \param record record of \a FileData_i in the file.
    \param FileData_i data read in the file at record \a record.
    \param FileData_f data read in the file at record \a record + 1.
    \param date date at which data is interpolated.
    \param InterpolatedData (output) data linearly interpolated between
    \a FileData_i and \a FileData_f.
  */
  template<class T>
  void BaseModel<T>::Interpolate(Date date_min_file, T Delta_t_file,
				 int record, const Data<T, 2>& FileData_i,
				 const Data<T, 2>& FileData_f, Date date,
				 Data<T, 2>& InterpolatedData)
  {
    int j, i;
    int Ny = FileData_i.GetLength(0);
    int Nx = FileData_i.GetLength(1);

    T time_distance = date.GetSecondsFrom(date_min_file);
    T weight_f = (time_distance - T(record) * Delta_t_file) / Delta_t_file;
    T weight_i = 1. - weight_f;

    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++)
	InterpolatedData(j, i) = FileData_i(j, i) * weight_i
	  + FileData_f(j, i) * weight_f;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_BASEMODEL_CXX
#endif
