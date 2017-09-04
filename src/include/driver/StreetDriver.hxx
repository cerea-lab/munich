// Copyright (C) 2016, CEREA - ENPC - EDF R&D
// Author(s): Youngseob Kim
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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

#ifndef POLYPHEMUS_FILE_DRIVER_STREETDRIVER_HXX

#include <map>
#include <vector>
#include <string>
#include "AtmoData.hxx"

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////
  // STREETDRIVER //
  //////////////////


  /*! \brief This class provides a driver for a street-network model that describes 
    the the atmospheric concentrations of pollutants in an urban street network.
  */
  /*! The driver is responsible for the model initialization, for the loop
    over all time steps, and for the calls to the output saver. 
    Its reference floating-point precision is 'T'. 
    The model is an instance of 'ClassModel' and the output saver is an
    instance of 'ClassOutputSaver'.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  class StreetDriver
  {

  protected:

    /*** Main components ***/

    //! Underlying model.
    ClassModel Model;
    //! Output saver.
    ClassOutputSaver OutputSaver;

    /*** Configuration ***/

    //! Configuration stream.
    ConfigStream config;

    //! Display options.
    map<string, bool> option_display;

  public:

    /*** Constructor and destructor ***/

    StreetDriver(string config_file);
    virtual ~StreetDriver();

    /*** Other methods ***/

    virtual void Run();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_STREETDRIVER_HXX
#endif
