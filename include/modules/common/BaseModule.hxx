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


#ifndef POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULE_HXX


namespace Polyphemus
{


  ////////////////
  // BASEMODULE //
  ////////////////


  //! This class should be the mother class of every module.
  /*! It defines the base interface of a module.
   */
  class BaseModule
  {

  public:

    template<class T>
    void Init(BaseModel<T>& Model);

    template<class T>
    void Forward(BaseModel<T>& Model);
    template<class T>
    void Forward_aer(BaseModel<T>& Model);

    template<class T>
    void Backward(BaseModel<T>& Model);

    bool IsRequired(string field);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULE_HXX
#endif
