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
  void SaverUnitStreetDryDeposition<T, ClassModel>::Init(ConfigStream& config_stream,
                                                         ClassModel& Model)
  {
    BaseSaverUnit<T, ClassModel>::Init(config_stream, Model);

    // Species to be saved.
    this->Ns = int(this->species_list.size());
    for (int s = 0; s < this->Ns; s++)
      this->species_index.push_back(Model.GetSpeciesIndex
				    (this->species_list[s]));

    Nstreet = Model.GetNStreet();
    StreetDryDepositionFlux_.Resize(this->Ns, Nstreet);
    WallDryDepositionFlux_.Resize(this->Ns, Nstreet);
    StreetDryDepositionRate_.Resize(this->Ns, Nstreet);
    WallDryDepositionRate_.Resize(this->Ns, Nstreet);

    // Output filename.
    string filename = config_stream.GetValue("Output_file");
    // Output filenames for all species.
    output_file.resize(this->Ns);
    output_file_wall.resize(this->Ns);
    output_file_rate.resize(this->Ns);
    output_file_rate_wall.resize(this->Ns);
    for (unsigned int i = 0; i < this->species_list.size(); i++)
      {
        output_file[i] = find_replace(filename, "&f", "FluxStreet_" + this->species_list[i]);
        output_file_wall[i] = find_replace(filename, "&f", "FluxWall_" + this->species_list[i]);
        output_file_rate[i] = find_replace(filename, "&f", "RateStreet_" + this->species_list[i]);
        output_file_rate_wall[i] = find_replace(filename, "&f", "RateWall_" + this->species_list[i]);
      }

    // Empties output files.
    for (unsigned int s = 0; s < this->species_list.size(); s++)
      {
        ofstream tmp_stream(output_file[s].c_str());
        ofstream tmp_stream_wall(output_file_wall[s].c_str());
        ofstream tmp_stream_rate(output_file_rate[s].c_str());
        ofstream tmp_stream_rate_wall(output_file_rate_wall[s].c_str());
      }

    // Initializing the buffer used to compute averaged concentrations.
    if (this->averaged)
      {
        int s, st;
        // StreetDryDepositionFlux_.Resize(this->Ns, Nstreet);        
        for (s = 0; s < this->Ns; s++)
          for (st = 0; st < Nstreet; st++)
            {
              StreetDryDepositionFlux_(s, st) = 0.5
                * Model.GetStreetDryDepositionFlux()(this->species_index[s], st);
              WallDryDepositionFlux_(s, st) = 0.5
                * Model.GetWallDryDepositionFlux()(this->species_index[s], st);
              StreetDryDepositionRate_(s, st) = 0.5
                * Model.GetStreetDryDepositionRate()(this->species_index[s], st);
              WallDryDepositionRate_(s, st) = 0.5
                * Model.GetWallDryDepositionRate()(this->species_index[s], st);
            }
      }

    // if (this->initial_concentration && !this->averaged)
    //   this->Save(Model);
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
                {
                  StreetDryDepositionFlux_(s, st) += 0.5
                    * Model.GetStreetDryDepositionFlux()
                    (this->species_index[s], st);
                  WallDryDepositionFlux_(s, st) += 0.5
                    * Model.GetWallDryDepositionFlux()
                    (this->species_index[s], st);
                  StreetDryDepositionRate_(s, st) += 0.5
                    * Model.GetStreetDryDepositionRate()
                    (this->species_index[s], st);
                  WallDryDepositionRate_(s, st) += 0.5
                    * Model.GetWallDryDepositionRate()
                    (this->species_index[s], st);
                }

            StreetDryDepositionFlux_.GetArray() /= T(this->interval_length);
            WallDryDepositionFlux_.GetArray() /= T(this->interval_length);
            StreetDryDepositionRate_.GetArray() /= T(this->interval_length);
            WallDryDepositionRate_.GetArray() /= T(this->interval_length);
            
            if (Model.GetCurrentDate() >= this->date_beg
        	&& Model.GetCurrentDate() <= this->date_end)
              for (s = 0; s < this->Ns; s++)
        	{
        	  Data<T, 1>
        	    StreetDryDepositionFlux_tmp(&StreetDryDepositionFlux_(s, 0), 
                                                shape(Nstreet));
        	  FormatBinary<float>().Append(StreetDryDepositionFlux_tmp,
        				       output_file[s]);

        	  Data<T, 1>
        	    WallDryDepositionFlux_tmp(&WallDryDepositionFlux_(s, 0), 
                                                shape(Nstreet));
        	  FormatBinary<float>().Append(WallDryDepositionFlux_tmp,
        				       output_file_wall[s]);

        	  Data<T, 1>
        	    StreetDryDepositionRate_tmp(&StreetDryDepositionRate_(s, 0), 
                                                shape(Nstreet));
        	  FormatBinary<float>().Append(StreetDryDepositionRate_tmp,
        				       output_file_rate[s]);

        	  Data<T, 1>
        	    WallDryDepositionRate_tmp(&WallDryDepositionRate_(s, 0), 
                                                shape(Nstreet));
        	  FormatBinary<float>().Append(WallDryDepositionRate_tmp,
        				       output_file_rate_wall[s]);
        	}


            for (s = 0; s < this->Ns; s++)
              for (st = 0; st < Nstreet; st++)
                {
                  StreetDryDepositionFlux_(s, st) = 0.5
                    * Model.GetStreetDryDepositionFlux()(this->species_index[s], st);
                  WallDryDepositionFlux_(s, st) = 0.5
                    * Model.GetWallDryDepositionFlux()(this->species_index[s], st);
                  StreetDryDepositionRate_(s, st) = 0.5
                    * Model.GetStreetDryDepositionRate()(this->species_index[s], st);
                  WallDryDepositionRate_(s, st) = 0.5
                    * Model.GetWallDryDepositionRate()(this->species_index[s], st);
                }
            this->counter = 0;
          }
        else
          {
            int s, st;
            for (s = 0; s < this->Ns; s++)
              for (st = 0; st < Nstreet; st++)
                {
                  StreetDryDepositionFlux_(s, st) +=
                    Model.GetStreetDryDepositionFlux()(this->species_index[s], st);
                  WallDryDepositionFlux_(s, st) +=
                    Model.GetWallDryDepositionFlux()(this->species_index[s], st);
                  StreetDryDepositionRate_(s, st) +=
                    Model.GetStreetDryDepositionRate()(this->species_index[s], st);
                  WallDryDepositionRate_(s, st) +=
                    Model.GetWallDryDepositionRate()(this->species_index[s], st);
                }
          }
      }
    else if (this->counter % this->interval_length == 0
             && Model.GetCurrentDate() >= this->date_beg
             && Model.GetCurrentDate() <= this->date_end)
      {
        // Instantaneous concentrations.
        for (int s = 0; s < this->Ns; s++)
          {
            Data<T, 1> StreetDryDepositionFlux_tmp(Nstreet);
            for(int st = 0; st < Nstreet; st++)
              StreetDryDepositionFlux_tmp(st) = Model.GetStreetDryDepositionFlux()
                (this->species_index[s], st);
            FormatBinary<float>().Append(StreetDryDepositionFlux_tmp, output_file[s]);

            Data<T, 1> WallDryDepositionFlux_tmp(Nstreet);
            for(int st = 0; st < Nstreet; st++)
              WallDryDepositionFlux_tmp(st) = Model.GetWallDryDepositionFlux()
                (this->species_index[s], st);
            FormatBinary<float>().Append(WallDryDepositionFlux_tmp, output_file_wall[s]);

            Data<T, 1> StreetDryDepositionRate_tmp(Nstreet);
            for(int st = 0; st < Nstreet; st++)
              StreetDryDepositionRate_tmp(st) = Model.GetStreetDryDepositionRate()
                (this->species_index[s], st);
            FormatBinary<float>().Append(StreetDryDepositionRate_tmp, output_file_rate[s]);

            Data<T, 1> WallDryDepositionRate_tmp(Nstreet);
            for(int st = 0; st < Nstreet; st++)
              WallDryDepositionRate_tmp(st) = Model.GetWallDryDepositionRate()
                (this->species_index[s], st);
            FormatBinary<float>().Append(WallDryDepositionRate_tmp, output_file_rate_wall[s]);
          }
      }
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREETDRYDEPOSITION_CXX
#endif
