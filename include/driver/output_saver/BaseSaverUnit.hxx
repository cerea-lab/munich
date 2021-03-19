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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_BASESAVERUNIT_HXX

#include <vector>
#include <fstream>

#include "AtmoDataHeader.hxx"

namespace Polyphemus
{


  //////////////
  // INCLUDES //
  //////////////


  using namespace std;
  using namespace AtmoData;


  ///////////////////
  // BASESAVERUNIT //
  ///////////////////


  //! This class is the base class for saver units.
  template<class T, class ClassModel>
  class BaseSaverUnit
  {

  protected:

    //! Origin of the domain along x.
    T base_x_min;
    //! Step along x.
    T base_Delta_x;
    //! Number of points along x.
    int base_Nx;
    //! Origin of the domain along y.
    T base_y_min;
    //! Step along y.
    T base_Delta_y;
    //! Number of points along y.
    int base_Ny;
    //! Number of vertical layers.
    int base_Nz;

    //! Counter of calls to the saver.
    int counter;

    //! Species list.
    vector<string> species_list;
    //! Number of species.
    int Ns;
    //! Species indices in the underlying model.
    vector<int> species_index;

    //! Starting date for the saver.
    Date date_beg;
    //! End date for the saver.
    Date date_end;
    //! Number of steps between two savings.
    int interval_length;
    //! Intantaneous or averaged concentrations.
    bool averaged;
    //! Save initial concentrations?
    bool initial_concentration;

  public:

    /*** Constructor and destructor ***/

    BaseSaverUnit();
    virtual ~BaseSaverUnit();

    /*** Initializations ***/

    virtual void Init(ConfigStream& config, ClassModel& Model);
    virtual void InitStep(ClassModel& Model);

    /*** Save ***/

    virtual void Save(ClassModel& Model) = 0;

    /*** Other methods ***/

    virtual string GetGroup() const;
    virtual string GetType() const = 0;

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_BASESAVERUNIT_HXX
#endif
