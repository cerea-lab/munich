// Copyright (C) 2007, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULE_CXX


#include "BaseModule.hxx"


namespace Polyphemus
{


  //! Initialization of the module.
  /*! \note Empty method.
   */
  template<class T>
  void BaseModule::Init(BaseModel<T>& Model)
  {
    throw "\"BaseModule::Init(BaseModel<T>& Model)\" is undefined.";
  }


  //! Performs an integration over one time step.
  /*!
    \note Empty method.
  */
  template<class T>
  void BaseModule::Forward(BaseModel<T>& Model)
  {
    throw "\"BaseModule::Forward(BaseModel<T>& Model)\" is undefined.";
  }


  //! Performs an integration over one time step for aerosols.
  /*!
    \note Empty method.
  */
  template<class T>
  void BaseModule::Forward_aer(BaseModel<T>& Model)
  {
    throw "\"BaseModule::Forward_aer(BaseModel<T>& Model)\" is undefined.";
  }


  //! Performs one step of backward integration of adjoint model.
  /*!
    \note Empty method.
  */
  template<class T>
  void BaseModule::Backward(BaseModel<T>& Model)
  {
    throw "\"BaseModule::Backward(BaseModel<T>& Model)\" is undefined.";
  }


  //! Checks whether a given field is required by this module.
  /*! Checks whether a given field must be available in the underlying model
    for this module to perform properly.
    \param field the field name.
    \note Always return false, must be redefined in derived classes.
  */
  bool BaseModule::IsRequired(string field)
  {
    return false;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULE_CXX
#endif
