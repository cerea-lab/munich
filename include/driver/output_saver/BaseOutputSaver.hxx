// Copyright (C) 2005-2016, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_BASEOUTPUTSAVER_HXX


#include "SaverUnitStreet.cxx"
#include "SaverUnitStreet_aer.cxx"
#include "SaverUnitStreetDryDeposition.cxx"
#include "SaverUnitStreetDryDepositionVelocity.cxx"
#include "SaverUnitStreetSurfaceDepositedMass.cxx"
#include "SaverUnitStreetSurfaceDepositedMass_aer.cxx"
#include "SaverUnitStreetResuspensionRate_aer.cxx"
#include "SaverUnitStreetWashoffRate_aer.cxx"
#include "SaverUnitStreetSurfaceDryDepositionRate_aer.cxx"
// #include "SaverUnitStreetDryDeposition.cxx"
// #include "SaverUnitStreetWetDeposition.cxx"

namespace Polyphemus
{


  //////////////
  // INCLUDES //
  //////////////


#include <vector>
#include <fstream>
  using namespace std;

#include "AtmoData.hxx"
  using namespace AtmoData;


  /////////////////////
  // BASEOUTPUTSAVER //
  /////////////////////


  //! This class manages a set of output-saver units.
  template<class T, class ClassModel>
  class BaseOutputSaver
  {

  protected:

    //! List of output-saver units.
    vector<BaseSaverUnit<T, ClassModel>*> Units;

    //! Group of active saver units.
    string group;

  public:

    BaseOutputSaver();
    ~BaseOutputSaver();

    void SetGroup(string grp);
    string GetGroup() const;
    bool OnlyCoordinatesList() const;
    
    void Init(ClassModel& Model);
    void InitStep(ClassModel& Model);
    void Save(ClassModel& Model);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_BASEOUTPUTSAVER_HXX
#endif
