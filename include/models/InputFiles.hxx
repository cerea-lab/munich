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


#ifndef POLYPHEMUS_FILE_MODELS_INPUTFILES_HXX


#include <vector>
#include <map>


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////
  // INPUTFILES //
  ////////////////


  //! This class stores a list of input files associated with a group of
  //! fields.
  /*! A group gathers several fields. Each field is associated with a path to
    the file that stores its values. In addition, the group is associated with
    a starting date and a time step that describes the files content.
  */
  template<class T>
  class InputFiles
  {

  protected:

    //! Starting date of files.
    Date date_min_;
    //! Time step of files.
    T Delta_t_;

    //! List of files associated to the fields.
    map<string, string> files_;

  public:

    /*** Constructor ***/

    InputFiles();

    /*** Basic access ***/

    T GetDelta_t() const;
    Date GetDateMin() const;

    string& operator()(string name);

    int GetLength() const;

    map<string, string>& GetMap();
    map<string, string>::iterator Begin();
    map<string, string>::iterator End();

    /*** Other methods ***/

    void Empty();
    void Read(string config_file, string section);
    void ReadFiles(string config_file, string section);
    void ReadFields(string config_file, string section);
    void Expand();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_INPUTFILES_HXX
#endif
