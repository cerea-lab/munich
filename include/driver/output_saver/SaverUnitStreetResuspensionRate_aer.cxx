// Copyright (C) 2016, CEREA - ENPC - EDF R&D
// Author(s): Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREETRESUSPENSIONRATE_AER_CXX


#include "SaverUnitStreetResuspensionRate_aer.hxx"


namespace Polyphemus
{


  /////////////////////
  // SAVERUNITSTREET //
  /////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitStreetResuspensionRate_aer<T, ClassModel>
  ::SaverUnitStreetResuspensionRate_aer():
    BaseSaverUnit<T, ClassModel>::BaseSaverUnit()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitStreetResuspensionRate_aer<T, ClassModel>::~SaverUnitStreetResuspensionRate_aer()
  {
  }


  //! Type of saver.
  template<class T, class ClassModel>
  string SaverUnitStreetResuspensionRate_aer<T, ClassModel>::GetType()  const
  {
    return "street_dep_mass_aer";
  }


  //! First initialization.
  /*! Reads the configuration.
    \param config_stream configuration stream.
    \param Model model with the following interface:
    <ul>
    <li> GetSpeciesIndex(string)
    <li> GetX_min()
    <li> GetDelta_x()
    <li> GetNx()
    <li> GetY_min()
    <li> GetDelta_y()
    <li> GetNy()
    <li> GetNz()
    <li> GetConcentration()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitStreetResuspensionRate_aer<T, ClassModel>::Init(ConfigStream& config_stream,
						ClassModel& Model)
  {

    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Check that memory has been allocated for StreetConcentration_aer.
    if(!Model.HasField("StreetResuspensionRate_aer"))
      throw string("The model does not have the array for street concentrations of particulate species.");

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species and bins.
    string field, species, file;
    vector<string> vsplit, bounds;
    int first, last, j, k;
    int i, Ns_aer;

    Ns_aer = Model.GetNs_aer();
    Nstreet = Model.D3("StreetResuspensionRate_aer").GetLength(2);
    
    vector<string> list_aer = Model.GetSpeciesList_aer();
    vector<int> bins;
    pair<string, vector<int> > tmp;
    for (i = 0; i < Model.GetNbin_aer(); i++)
      bins.push_back(i);
    tmp.second = bins;
    for (i = 0; i < Ns_aer; i++)
      {
	tmp.first = list_aer[i];
	species_list_aer.push_back(tmp);
	// Adds a buffer to compute averaged concentrations.
	Concentration_aer_.push_back(vector<Data<T, 1> >());
	output_file.push_back(vector<string>());
	for (k = 0; k < int(bins.size()); k++)
	  {
	    file = find_replace(filename, "&f", tmp.first);
	    file = find_replace(file, "&n", to_str(bins[k]));
	    output_file[i].push_back(file);
	    Concentration_aer_[i].push_back(Data<T, 1>());
	  }
      }

    // Indices in the model.
    int base_s, base_b;

    int s, b, st;
    // Empties output files.
    for (s = 0; s < int(species_list_aer.size()); s++)
      for (b = 0; b < int(species_list_aer[s].second.size()); b++)
        ofstream tmp_stream(output_file[s][b].c_str());
    
    // Initializing the buffer used to compute averaged concentrations.
    if (this->averaged)
      for (s = 0; s < int(species_list_aer.size()); s++)
	for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	  {
	    base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
	    base_b = species_list_aer[s].second[b];
	    Concentration_aer_[s][b].Resize(Nstreet);
	    for (st = 0; st < Nstreet; st++)
	      Concentration_aer_[s][b](st) = 0.5
		* Model.D3("StreetResuspensionRate_aer")(base_s, base_b, st);
	      
	  }

    if (this->initial_concentration && !this->averaged)
      this->Save(Model);
  }

  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitStreetResuspensionRate_aer<T, ClassModel>::InitStep(ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::InitStep(Model);
  }


  //! Saves concentrations if needed.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetConcentration()
    <li> GetCurrentDate()
    <li> GetConcentration(int, T, T, T)
    <li> GetIntegratedConcentration(int, T, T, T, T, T, T)
    <li> ComputeConcentration()
    </ul>
  */
  template<class T, class ClassModel>
  void SaverUnitStreetResuspensionRate_aer<T, ClassModel>::Save(ClassModel& Model)
  {
    int s, b, k, j, i, st;
    // Indices in the model.
    int base_s, base_b;
    base_s = -999;
    
    if (this->averaged)
      if (this->counter % this->interval_length == 0)
	for (s = 0; s < int(species_list_aer.size()); s++)
	  for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	    {
	      base_b = species_list_aer[s].second[b];
	      
	      base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
	      for (st = 0; st < Nstreet; st++)
		Concentration_aer_[s][b](st) += 0.5
		  * Model.D3("StreetResuspensionRate_aer")(base_s, base_b, st);
		
	      Concentration_aer_[s][b].GetArray() /= T(this->interval_length);
	      
	      if (Model.GetCurrentDate() >= this->date_beg
		  && Model.GetCurrentDate() <= this->date_end)
		FormatBinary<float>().Append(Concentration_aer_[s][b],
					     output_file[s][b]);
	      
	      for (st = 0; st < Nstreet; st++)
		Concentration_aer_[s][b](st) = 0.5
		  * Model.D3("StreetResuspensionRate_aer")(base_s, base_b, st);
	
	      this->counter = 0;
	    }
    //
      else
	for (s = 0; s < int(species_list_aer.size()); s++)
	  for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	    {
	      base_b = species_list_aer[s].second[b];
	      base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
	      for (st = 0; st < Nstreet; st++)
		Concentration_aer_[s][b](st) +=
		  Model.D3("StreetResuspensionRate_aer")(base_s, base_b, st);
	
	    }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      {
	for (s = 0; s < int(species_list_aer.size()); s++)
	  for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	    {
	      base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
	      base_b = species_list_aer[s].second[b];
	      // for (st = 0; st < Nstreet; st++)
	      // 	{
	      Data<T, 1> Concentration_tmp(&Model.D3("StreetResuspensionRate_aer")
					   (base_s, base_b,0),
					   shape(Nstreet));
	      FormatBinary<float>().Append(Concentration_tmp,
					   output_file[s][b]);
	      //		}
	    }
      }
    
  } // namespace Polyphemus.
}

#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREETRESUSPENSIONRATE_AER_CXX
#endif
