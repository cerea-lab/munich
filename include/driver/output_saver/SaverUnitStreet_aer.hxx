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

//#ifdef POLYPHEMUS_WITH_AEROSOL_MODULE
#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREET_AER_HXX


#include "BaseSaverUnit.cxx"
// #include "SaverUnitStreet.cxx"
#include <vector>
#include <fstream>
#include "AtmoData.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  /////////////////////
  // SAVERUNITSTREET //
  /////////////////////


  /*! \brief This class is used to save concentrations over the whole domain
    of the underlying model. Chemical species and vertical levels may be
    selected. */
  template<class T, class ClassModel>
  class SaverUnitStreet_aer: public BaseSaverUnit<T, ClassModel>
  {

  protected:

    //! List of aerosol species (with their bins) to be saved.
    vector<pair<string, vector<int> > > species_list_aer;

    //! List of output files.
    vector<vector<string> > output_file;
    
    /*! Buffer used to compute averaged concentrations. It is first indexed by
      the species name and the bin. */
    vector<vector<Data<T, 2> > > Concentration_;
    
    //! Number of points to be saved.
    int Nstreet;  
    
    //! Type of saver (indices or coordinates).
    string type;    

  public:

    SaverUnitStreet_aer();
    virtual ~SaverUnitStreet_aer();
    
    virtual string GetType() const;

    virtual void Init(ConfigStream& config_stream, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);
    virtual void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_SAVERUNITSTREET_AER_HXX
#endif
//#endif
