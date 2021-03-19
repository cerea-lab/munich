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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREETDRYDEPOSITION_CXX


#include "SaverUnitStreetDryDeposition.hxx"


namespace Polyphemus
{


  //////////////////////////////////
  // SAVERUNITSTREETDRYDEPOSITION //
  //////////////////////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  SaverUnitStreetDryDeposition<T, ClassModel>
  ::SaverUnitStreetDryDeposition():
    BaseSaverUnit<T, ClassModel>::BaseSaverUnit()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  SaverUnitStreetDryDeposition<T, ClassModel>::~SaverUnitStreetDryDeposition()
  {
  }


  //! Type of saver.
  template<class T, class ClassModel>
  string SaverUnitStreetDryDeposition<T, ClassModel>::GetType()  const
  {
    return "street_dry_deposition";
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
  void SaverUnitStreetDryDeposition<T, ClassModel>
  ::Init(ConfigStream& config_stream,
	 ClassModel& Model)
  {

    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Check that memory has been allocated for StreetConcentration.
    if(!Model.HasField("StreetDryDepositionRate"))
      throw string("The model does not have the array for street dry deposition rates.");

    // Species to be saved.
    this->Ns = int(this->species_list.size());
    vector<string> model_species_list = Model.GetSpeciesList();
    vector<string> deposition_species_list =
      Model.GetSpeciesList("DepositionVelocity");
    int deposition_index;
    
    for (int s = 0; s < this->Ns; s++)
      {
	if (find(model_species_list.begin(), model_species_list.end(),
		 this->species_list[s]) == model_species_list.end())
	  throw string("Species \"") + this->species_list[s] + "\" unknown.";
	else if (find(deposition_species_list.begin(),
		      deposition_species_list.end(), this->species_list[s])
		 == deposition_species_list.end())
	  throw string("Species \"") + this->species_list[s]
	    + "\" is not dry deposited.";
	else
	  {
	    deposition_index = 0;
	    while (deposition_species_list[deposition_index] != this->species_list[s])
	      deposition_index++;
	    this->species_index.push_back(deposition_index);
	  }
      }

    //    Nstreet = Model.GetNStreet();
    Nstreet = Model.D2("StreetDryDepositionRate").GetLength(1);
    StreetDryDepositionRate_.Resize(this->Ns, Nstreet);

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species.
    output_file.resize(this->Ns);
    for (unsigned int i = 0; i < this->species_list.size(); i++)
      output_file[i] = find_replace(filename, "&f", this->species_list[i]);

    // Empties output files.
    for (unsigned int s = 0; s < this->species_list.size(); s++)
      ofstream tmp_stream(output_file[s].c_str());

    // Initializing the buffer used to compute averaged concentrations.
    if (this->averaged)
      {
        int s, st;
        StreetDryDepositionRate_.Resize(this->Ns, Nstreet);
        for (s = 0; s < this->Ns; s++)
          for (st = 0; st < Nstreet; st++)
            StreetDryDepositionRate_(s, st) = 0.5
              * Model.D2("StreetDryDepositionRate")(this->species_index[s], st);

      }
  }

  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void SaverUnitStreetDryDeposition<T, ClassModel>::InitStep(ClassModel& Model)
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
  void SaverUnitStreetDryDeposition<T, ClassModel>::Save(ClassModel& Model)
  {
    if (this->averaged)
      {
        if (this->counter % this->interval_length == 0)
          {
            int s, st;
            for (s = 0; s < this->Ns; s++)
	      for (st = 0; st < Nstreet; st++)
                StreetDryDepositionRate_(s, st) += 0.5
                  * Model.D2("StreetDryDepositionRate")(this->species_index[s], st);

            StreetDryDepositionRate_.GetArray() /= T(this->interval_length);
            
            if (Model.GetCurrentDate() >= this->date_beg
        	&& Model.GetCurrentDate() <= this->date_end)
              for (s = 0; s < this->Ns; s++)
        	{
        	  Data<T, 1>
        	    StreetDryDepositionRate_tmp(&StreetDryDepositionRate_(s, 0), shape(Nstreet));
        	  FormatBinary<float>().Append(StreetDryDepositionRate_tmp,
        				       output_file[s]);
        	}


            for (s = 0; s < this->Ns; s++)
              for (st = 0; st < Nstreet; st++)
                StreetDryDepositionRate_(s, st) = 0.5
                  * Model.D2("StreetDryDepositionRate")(this->species_index[s], st);

            this->counter = 0;
          }
        else
          {
            int s, st;
            for (s = 0; s < this->Ns; s++)
              for (st = 0; st < Nstreet; st++)
                StreetDryDepositionRate_(s, st) +=
                  Model.D2("StreetDryDepositionRate")(this->species_index[s], st);
          }
      }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      {
        // Instantaneous deposition rate.
        for (int s = 0; s < this->Ns; s++)
          {
            Data<T, 1> DepositionRate_street(Nstreet);
            for(int st = 0; st < Nstreet; st++)
              DepositionRate_street(st) = Model.D2("StreetDryDepositionRate")(this->species_index[s], st);
            FormatBinary<float>().Append(DepositionRate_street, output_file[s]);
          }
      }

  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREETDEPOSITIONRATE_CXX
#endif
