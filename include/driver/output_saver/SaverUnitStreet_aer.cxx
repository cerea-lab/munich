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
   
    // Checks whether the number is computed in the model.
    // Use the virtual function "HasNumberConcentration_aer" of BaseModel  
    // Total number of streets.
    Nstreet = Model.GetNStreet();

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
		      * Model.GetStreetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition, st);
	      }
	    else if (model_with_number)
	      {
		base_b = species_list_aer[s].second[b];
		Concentration_[s][b].Resize(Model.GetNcomposition_aer(), Nstreet);
		 for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) = 0.5
			* Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition, st);
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
			* Model.GetStreetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition, st);
		}
	      else if (model_with_number)		
		{
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition, st) += 0.5
			* Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition, st);
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
			* Model.GetStreetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition, st);
		}
	      else if (model_with_number)
		{
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) = 0.5
			* Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition, st);
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
			Model.GetStreetConcentration_aer()(base_s, base_b*Model.GetNcomposition_aer()+icomposition, st);
		}
	      else if (model_with_number)
		{
		  for(icomposition=0; icomposition< Model.GetNcomposition_aer();icomposition++)
		    for (int st = 0; st < Nstreet; st++)
		      Concentration_[s][b](icomposition,st) +=
			Model.GetStreetNumberConcentration()(base_b*Model.GetNcomposition_aer()+icomposition, st);
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
