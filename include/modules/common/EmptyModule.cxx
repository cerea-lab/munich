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


#ifndef POLYPHEMUS_FILE_MODULES_COMMON_EMPTYMODULE_CXX


#include "EmptyModule.hxx"


namespace Polyphemus
{


  //! Initialization of the module.
  /*!
    \note Empty method.
  */
  template<class ClassModel>
  void EmptyModule::Init(ClassModel& Model)
  {
  }


  //! Performs an integration over one time step.
  /*!
    \note Empty method.
  */
  template<class ClassModel>
  void EmptyModule::Forward(ClassModel& Model)
  {
  }


  //! Performs an integration over one time step for aerosols.
  /*!
    \note Empty method.
  */
  template<class ClassModel>
  void EmptyModule::Forward_aer(ClassModel& Model)
  {
  }


  //! Performs one step of backward integration of adjoint model.
  /*!
    \note Empty method.
  */
  template<class ClassModel>
  void EmptyModule::Backward(ClassModel& Model)
  {
  }


  //! Checks whether a given field is required by this module.
  /*! Checks whether a given field must be available in the underlying model
    for this module to perform properly.
    \param[in] field the field name.
    \note Always returns false.
  */
  bool EmptyModule::IsRequired(string field)
  {
    return false;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_COMMON_EMPTYMODULE_CXX
#endif
