// Copyright (C) 2009 INRIA
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


#ifndef POLYPHEMUS_FILE_MODULES_COMMON_EMPTYMODULE_HXX


namespace Polyphemus
{


  /////////////////
  // EMPTYMODULE //
  /////////////////


  //! This class defines an empty module.
  /*! This module performs nothing at all!
   */
  class EmptyModule
  {

  public:

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void Forward(ClassModel& Model);
    template<class ClassModel>
    void Forward_aer(ClassModel& Model);

    template<class ClassModel>
    void Backward(ClassModel& Model);

    bool IsRequired(string field);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_COMMON_EMPTYMODULE_HXX
#endif
