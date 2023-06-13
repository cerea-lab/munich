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


#ifndef POLYPHEMUS_FILE_DRIVER_STREETDRIVER_CXX


#include "StreetDriver.hxx"

namespace Polyphemus
{

  //! Main constructor.
  /*! Builds the driver and reads option keys in the configuration file.
    \param config_file configuration file.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  StreetDriver<T, ClassModel, ClassOutputSaver>
  ::StreetDriver(string config_file):
    Model(config_file), config(config_file)
  {

    /*** Display options ***/

    config.SetSection("[display]");
    // Should iterations be displayed on screen?
    config.PeekValue("Show_iterations", option_display["show_iterations"]);
    // Should current date be displayed on screen?
    config.PeekValue("Show_date", option_display["show_date"]);
  }

  //! Destructor.
  template<class T, class ClassModel, class ClassOutputSaver>
  StreetDriver<T, ClassModel, ClassOutputSaver>::~StreetDriver()
  {
  }

  //! Performs the simulation.
  /*! Initializes the model and the output saver, and then performs the loop
    over all meteorological conditions with calls to the model (over all time
    steps) and to the output saver.
  */
  template<class T, class ClassModel, class ClassOutputSaver>
  void StreetDriver<T, ClassModel, ClassOutputSaver>::Run()
  {

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    int argc; 
    char **argv; 

    // YK    MPI::Init();
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    string line;

    /*** Initializations ***/

    //read configuration options in configuration files
    Model.ReadConfiguration();
    //Read number of streets and intersections
    Model.Init();
    //Initialize input data
    Model.InitAllData();
    
    OutputSaver.Init(Model);
    Model.InitOutputSaver();
    
    for (int i = 0; i < Model.GetNt(); i++)
      {
	if (option_display["show_iterations"] && (rank == 0))
	  cout << "Performing iteration #" << i << endl;

	if (option_display["show_date"] && (rank == 0))
	  cout << "Current date: "
	       << Model.GetCurrentDate().GetDate("%y-%m-%d %h:%i") << endl;
	//update input data
	Model.InitStep();
	if (rank == 0)
	  OutputSaver.InitStep(Model);

	Model.Forward();

	if (rank == 0)
          {
	    OutputSaver.Save(Model);
	    Model.OutputSaver();
	  }
      }
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // YK    MPI::Finalize();
    MPI_Finalize();    
#endif
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_STREETDRIVER_CXX
#endif
