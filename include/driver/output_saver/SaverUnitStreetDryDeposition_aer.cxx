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

#ifdef POLYPHEMUS_WITH_AEROSOL_MODULE
#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREET_AER_CXX


#include "SaverUnitStreet_aer.hxx"


namespace Polyphemus
{


  /////////////////////
  // SAVERUNITSTREET //
  /////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitStreet_aer<T, ClassModel>
  ::SaverUnitStreet_aer():
    BaseSaverUnit<T, ClassModel>::BaseSaverUnit()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitStreet_aer<T, ClassModel>::~SaverUnitStreet_aer()
  {
  }


  //! Type of saver.
  template<class T, class ClassModel>
  string SaverUnitStreet_aer<T, ClassModel>::GetType()  const
  {
    return type;
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
  void SaverUnitStreet_aer<T, ClassModel>::Init(ConfigStream& config_stream,
                                            ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // // Check that memory has been allocated for StreetConcentration.
    // if(!Model.HasField("StreetConcentration_aer"))
    //   throw string("The model does not have the array for aerosol street concentrations.");
    
    // Total number of streets.
    Nstreet = Model.GetNStreet();

    // Checks whether the number is computed in the model.
    // Use the virtual function "HasNumberConcentration_aer" of BaseModel  
    bool model_with_number = Model.HasNumberConcentration_aer();
    //cout<<"saver 0"<<endl;

    // Output filename.
    string filename = config_stream.GetValue("Output_file");  
    // Output filenames for all species and bins.
    string field, species, file;
    vector<string> vsplit, bounds;
    int first, last, j, k;
    
    int Ns_aer = Model.GetNs_aer();
    int Ns_aer_all = Model.GetNs_aer();
    vector<string> list_aer = Model.GetSpeciesList_aer();
    
    if (model_with_number)
      {
	list_aer.push_back("Number");
	Ns_aer++;
      }    

    // if (this->species_list[0] == "all")//all aerosol species
    //   {
    int i;
    vector<int> bins;
    pair<string, vector<int> > tmp;
    // KS Check that Model.GetNbin_aer(GetNsize_section_aer) output something between 0 and Nsize_section
    for (i = 0; i < Model.GetNsize_section_aer(); i++)
      bins.push_back(i);

    tmp.second = bins;
    for (i = 0; i < Ns_aer; i++)
      {
	tmp.first = list_aer[i];
	species_list_aer.push_back(tmp);
	// Adds a buffer to compute averaged concentrations.
	Concentration_.push_back(vector<Data<T, 2> >());
	output_file.push_back(vector<string>());
	for (k = 0; k < int(bins.size()); k++)
	  {
	    file = find_replace(filename, "&f", tmp.first);
	    file = find_replace(file, "&n", to_str(bins[k]));	
	    output_file[i].push_back(file);
	    Concentration_[i].push_back(Data<T, 2>()); 
	  }
      }  
	//}
      
    //else
      // {
      // 	for (unsigned int i = 0; i < this->species_list.size(); i++)
      // 	  {
      // 	    field = this->species_list[i];
      // 	    vsplit = split(field, "{}");
      // 	    if (field[0] == '{')
      // 	      throw string("Species \"") + field + string("\" is badly ")
      // 		+ "formatted: it cannot be parsed by the output saver.";
      // 	    if (vsplit.size() == 1)
      // 	      throw string("Species \"") + field + string("\" is badly ")
      // 		+ "formatted: bins are needed by the output saver.";
      // 	    if (vsplit.size() > 2)
      // 	      throw string("Species \"") + field + string("\" is badly ")
      // 		+ "formatted: it cannot be parsed by the output saver.";

      // 	    // Searches for species index in 'species_list_aer' and
      // 	    // 'output_file'.  If the species is not found, it is added in
      // 	    // 'species_list_aer' and 'output_file'.
      // 	    species = split(field, "_")[0];
      // 	    j = 0;
		  
      // 	    while (j < int(species_list_aer.size())
      // 		   && species_list_aer[j].first != species)
      // 	      j++;
      // 	    if (j == int(species_list_aer.size()))
      // 	      {
      // 		species_list_aer.push_back(pair<string, vector<int> >
      // 					   (species, vector<int>()));
      // 		output_file.push_back(vector<string>());
      // 		Concentration_.push_back(vector<Data<T, 2> >()); 
      // 	      }
      // 	    bounds = split(vsplit[1], "-");
      // 	    // First bound.
      // 	    first = convert<int>(bounds[0]);
      // 	    // Last bound.
      // 	    if (bounds.size() != 1)
      // 	      last = convert<int>(bounds[1]);
      // 	    else
      // 	      last = first;
      // 	    for (k = first; k < last + 1; k++)
      // 	      {
      // 		// Adds the bin (associated to species #j).
      // 		species_list_aer[j].second.push_back(k);
      // 		// Adds a buffer to compute averaged concentrations.
      // 		Concentration_[j].push_back(Data<T, 2>()); 

      // 		// Field base name is replaced in the generic file name.
      // 		file = find_replace(filename, "&f",
      // 				    species);
     
      // 		// Field number is also replaced.
      // 		file = find_replace(file, "&n", to_str(k));
      // 		// Adds the corresponding output file.
      // 		output_file[j].push_back(file);
      // 	      }
      // 	  }
	
      // 	if (model_with_number)
      // 	  {
      // 	    output_file.push_back(vector<string>());
      // 	    Concentration_.push_back(vector<Data<T, 2> >()); 
		
      // 	    vector<int> bins;
      // 	    pair<string, vector<int> > tmp;
      // 	    for (int i = 0; i < Model.GetNsize_section_aer(); i++)
      // 	      bins.push_back(i);
      // 	    tmp.second = bins;
      // 	    tmp.first = "Number";
      // 	    species_list_aer.push_back(tmp);
			
      // 	    int rear = output_file.size()-1;
      // 	    for (k = 0; k < int(bins.size()); k++)
      // 	      {
      // 		file = find_replace(filename, "&f", "Number");
      // 		file = find_replace(file, "&n", to_str(bins[k]));
      // 		output_file[rear].push_back(file);
      // 		Concentration_[rear].push_back(Data<T, 2>());
      // 	      }
      // 	  }
      // }
    //cout<<"saver 1"<<endl;
    // Indices in the model.
    int base_s, base_b;

    int s, b, icomposition;
    // Empties output files.
    for (s = 0; s < int(species_list_aer.size()); s++)
      for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	ofstream tmp_stream(output_file[s][b].c_str());
    
    if (this->averaged)
      for (s = 0; s < int(species_list_aer.size()); s++)
	for (b = 0; b < int(species_list_aer[s].second.size()); b++)	  
	  {   
	    if (species_list_aer[s].first != "Number")
	      {
		base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
		base_b = species_list_aer[s].second[b];
		Concentration_[s][b].Resize(Model.GetNcomposition_aer(), Nstreet);
	        for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		  for (int st = 0; st < Nstreet; st++)
		    Concentration_[s][b](icomposition, st) = 0.5
		      * Model.GetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition);
	      }
	    else if (model_with_number)
	      {
		base_b = species_list_aer[s].second[b];
		Concentration_[s][b].Resize(Model.GetNcomposition_aer(), Nstreet);
		 for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) = 0.5
			  * Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition);
	      }
	  }
    //cout<<"saver 2"<<endl;
    if (this->initial_concentration && !this->averaged)
      Save(Model);
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitStreet_aer<T, ClassModel>::InitStep(ClassModel& Model)
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
  void SaverUnitStreet_aer<T, ClassModel>::Save(ClassModel& Model)
  {
    int s, b, k, j, i,icomposition;
    // Indices in the model.
    int base_s, base_b;
    bool model_with_number = Model.HasNumberConcentration_aer();
    
    // Total number of streets.
    Nstreet = Model.GetNStreet();
    if(this->counter % this->interval_length == 0)
      cout<<"Saving aerosol"<<this->counter<<"/"<<this->interval_length<<endl;
    else
      cout<<"Skip saving aerosol"<<this->counter<<"/"<<this->interval_length<<endl;
//       //save vertical winds
//       string File_out = "vertical_wind.bin";
//       FormatBinary<float>().Append(Model.GetVertical_Wind(), File_out);
    
    if (this->averaged)
      if (this->counter % this->interval_length == 0)
	for (s = 0; s < int(species_list_aer.size()); s++)
	  for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	    {
	      base_b = species_list_aer[s].second[b];
		  
	      if (species_list_aer[s].first != "Number")
		{
		  base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
		 for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) += 0.5
			* Model.GetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition);
		}
	      else if (model_with_number)		
		{
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition, st) += 0.5
			* Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition);
		}

	      Concentration_[s][b].GetArray() /= T(this->interval_length);

	      if (Model.GetCurrentDate() >= this->date_beg
		  && Model.GetCurrentDate() <= this->date_end)
	      {
		//cout<<output_file[s][b]<<endl;
		FormatBinary<float>().Append(Concentration_[s][b],
					     output_file[s][b]);
	      }
	      if (species_list_aer[s].first != "Number")
		{
		  base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) = 0.5
			* Model.GetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition);
		}
	      else if (model_with_number)
		{
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) = 0.5
			* Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition);
		}

	      this->counter = 0;
	    }
      else
	for (s = 0; s < int(species_list_aer.size()); s++)
	  for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	    {
	      base_b = species_list_aer[s].second[b];
	      if (species_list_aer[s].first != "Number")
		{
		  base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) +=
			Model.GetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition);
		}
	      else if (model_with_number)
		{
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) +=
			Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition);
		}
	    }
	
    else if (this->counter % this->interval_length == 0
	     &&	Model.GetCurrentDate() >= this->date_beg
	     && Model.GetCurrentDate() <= this->date_end)
      // Instantaneous concentrations.
      for (s = 0; s < int(species_list_aer.size()); s++)
	for (b = 0; b < int(species_list_aer[s].second.size()); b++)
	  {
	    base_b = species_list_aer[s].second[b];
	    for (icomposition=0;icomposition < Model.GetNcomposition_aer(); icomposition++)
	      if (species_list_aer[s].first != "Number")
		{
		  base_s = Model.GetSpeciesIndex_aer(species_list_aer[s].first);

		  Data<T, 1> Concentration_tmp(&Model.GetStreetConcentration_aer()
					      (base_s, base_b*Model.GetNcomposition_aer()+icomposition,0),
					      shape(Nstreet));
		  FormatBinary<float>().Append(Concentration_tmp,
					      output_file[s][b]);
// 		  cout<<output_file[s][b]<<endl;
		}
	      else if (model_with_number)
		{
		  Data<T, 1> NumberConcentration_tmp(&Model.GetStreetNumberConcentration()
						    (base_b*Model.GetNcomposition_aer()+icomposition,0),
						    shape(Nstreet));
		  FormatBinary<float>().Append(NumberConcentration_tmp,
					      output_file[s][b]);
// 		  cout<<output_file[s][b]<<endl;
		}

	  }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREET_AER_CXX
#endif
#endif
