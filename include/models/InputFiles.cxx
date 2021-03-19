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


#ifndef POLYPHEMUS_FILE_MODELS_INPUTFILES_CXX


#include "InputFiles.hxx"
#include "AtmoData.hxx"


namespace Polyphemus
{


  /////////////////
  // CONSTRUCTOR //


  //! Default constructor.
  template<class T>
  InputFiles<T>::InputFiles()
  {
    this->Empty();
  }


  // CONSTRUCTOR //
  /////////////////


  ////////////
  // ACCESS //


  //! Returns the files starting-date.
  /*! \return The files starting-date.
   */
  template<class T>
  Date InputFiles<T>::GetDateMin() const
  {
    return date_min_;
  }


  //! Returns the files time-step.
  /*! \return The files time-step.
   */
  template<class T>
  T InputFiles<T>::GetDelta_t() const
  {
    return Delta_t_;
  }


  //! Returns the path to the file associated with a given field.
  /*!
    \param name field name.
    \return The path to the file associated with the field 'name'.
  */
  template<class T>
  string& InputFiles<T>::operator()(string name)
  {
    if (files_[name] == "")
      cout << "Warning: undefined variable " << name << ".\n"; 
    return files_[name];
  }


  //! Returns the number of files.
  /*! \return The number of files.
   */
  template<class T>
  int InputFiles<T>::GetLength() const
  {
    return files_.size();
  }


  //! Returns the map between field names and file paths.
  /*! \return The map between field names and file paths.
   */
  template<class T>
  map<string, string>& InputFiles<T>::GetMap()
  {
    return files_;
  }


  /*! \brief Returns an iterator to the first element of the map between field
    names and file paths. */
  /*! \return An iterator to the first element of the map between field names
    and file paths.
  */
  template<class T>
  map<string, string>::iterator InputFiles<T>::Begin()
  {
    return files_.begin();
  }


  /*! \brief Returns the "end iterator" of the map between field names and
    file paths. */
  /*! \return The "end iterator" of the map between field names and file
    paths.
  */
  template<class T>
  map<string, string>::iterator InputFiles<T>::End()
  {
    return files_.end();
  }


  // ACCESS //
  ////////////


  ///////////////////
  // OTHER METHODS //


  //! Clears the object.
  /*! The starting date is arbitrarily set to 0000-01-01, and the time step to
    0.
  */
  template<class T>
  void InputFiles<T>::Empty()
  {
    date_min_.SetDate(101);
    Delta_t_ = 0;
    files_.clear();
  }


  /*! \brief Reads the starting date, the time step and all file paths from a
    configuration file. */
  /*! The configuration file must contain a section called [\a section] in
    which the starting date is the field "Date_min" and the time step (in
    seconds) is the field "Delta_t". The starting date should be in any format
    handled by the library Talos. At the end (this is compulsory) of the
    section, field names are put after "Fields". If no field name is provided,
    "Fields" must be followed by "-", "--" or "---". Next, a corresponding
    file path is provided in the (configuration-file) field "Filename". This
    is a generic path name in which any occurrence of "&f" is replaced with
    the field name. Finally, an additional list of field names and file paths
    may be provided. Example:

    \verbatim
    # This is a comment line.
    Date_min: 2006-06-13_12h30
    Delta_t: 600.

    Fields: Pressure Temperature
    # Paths: /home/user/data/meteo/Pressure.bin
    #        and /home/user/data/meteo/Temperature.bin
    Filename: /home/user/data/meteo/&f.bin

    # Other Fields.
    Rain   /home/user/data/new_meteo/Rain.bin
    Albedo /home/user/data/ground/RawAlbedo.bin \endverbatim

    \param config_file configuration file.
    \param section section in the configuration file.
  */
  template<class T>
  void InputFiles<T>::Read(string config_file, string section)
  {
    ConfigStream config_stream(config_file);
    config_stream.SetSection(string("[") + section + string("]"));

    date_min_ = config_stream.PeekValue("Date_min");
    config_stream.PeekValue("Delta_t", "> 0", Delta_t_);
    config_stream.Find("Fields");
    vector<string> fields;
    string fields_string = config_stream.GetLine();
    // If there is any field.
    if (trim(fields_string) != "---" && trim(fields_string) != "--"
        && trim(fields_string) != "-")
      fields = split(fields_string);

    string generic_filename = config_stream.GetValue("Filename");

    for (unsigned int i = 0; i < fields.size(); i++)
      files_[fields[i]] = generic_filename;

    // Other fields (not specified in the configuration-file field 'Fields').
    string field;
    while (!config_stream.IsEmpty())
      {
        field = config_stream.GetElement();
        files_[field] = config_stream.GetElement();
      }

    Expand();
  }


  //! Reads all file paths from a configuration file.
  /*! The configuration file must contain a section called [\a section] in
    which field names are put after "Fields". If no field name is provided,
    "Fields" must be followed by "-", "--" or "---". Next, a corresponding
    file path is provided in the (configuration-file) field "Filename". This
    is a generic path name in which any occurrence of "&f" is replaced with
    the field name. Finally, an additional list of field names and file paths
    may be provided. Example:

    \verbatim
    # This is a comment line.
    Fields: Pressure Temperature
    # Paths: /home/user/data/meteo/Pressure.bin
    #        and /home/user/data/meteo/Temperature.bin
    Filename: /home/user/data/meteo/&f.bin

    # Other Fields.
    Rain   /home/user/data/new_meteo/Rain.bin
    Albedo /home/user/data/ground/RawAlbedo.bin \endverbatim

    \param config_file configuration file.
    \param section section in the configuration file.
  */
  template<class T>
  void InputFiles<T>::ReadFiles(string config_file, string section)
  {
    ConfigStream config_stream(config_file);
    config_stream.SetSection(string("[") + section + string("]"));

    config_stream.Find("Fields");
    vector<string> fields;
    string fields_string = config_stream.GetLine();
    // If there is any field.
    if (trim(fields_string) != "---" && trim(fields_string) != "--"
        && trim(fields_string) != "-")
      fields = split(fields_string);

    string generic_filename = config_stream.GetValue("Filename");

    for (unsigned int i = 0; i < fields.size(); i++)
      files_[fields[i]] = generic_filename;

    // Other fields (not specified in the configuration-file field 'Fields').
    string field;
    while (!config_stream.IsEmpty())
      {
        field = config_stream.GetElement();
        files_[field] = config_stream.GetElement();
      }

    Expand();
  }


  //! Reads fields names from a configuration file.
  /*! The configuration file must contain a section called [\a section] in
    which field names are put after "Fields". If no field name is provided,
    "Fields" must be followed by "-", "--" or "---". Example:

    \verbatim
    # This is a comment line.
    Fields: Pressure Temperature Rain Albedo
    \endverbatim

    \param config_file configuration file.
    \param section section in the configuration file.
  */
  template<class T>
  void InputFiles<T>::ReadFields(string config_file, string section)
  {
    ConfigStream config_stream(config_file);
    config_stream.SetSection(string("[") + section + string("]"));

    config_stream.Find("Fields");
    vector<string> fields;
    string fields_string = config_stream.GetLine();
    if (trim(fields_string) != "---" && trim(fields_string) != "--"
        && trim(fields_string) != "-")
      fields = split(fields_string);
    for (unsigned i = 0; i < fields.size(); i++)
      files_[fields[i]] = "---";

    Expand();
  }


  //! Expands numbers in field names and file names.
  template<class T>
  void InputFiles<T>::Expand()
  {
    map<string, string> old_files = files_;
    files_.clear();

    map<string, string>::iterator iter;
    string field, file;
    vector<string> vsplit, bounds;
    int first, last;
    // Loop over all fields and associated files.
    for (iter = old_files.begin(); iter != old_files.end(); iter++)
      {
        // Field name.
        field = iter->first;
        // File name.
        file = iter->second;

        // Splits the field name and its numbers (enclosed in {}).
        vsplit = split(field, "{}");
        if (field[0] == '{')
          // Only numbers (no field base name).
          {
            bounds = split(vsplit[0], "-");
            // First bound.
            first = convert<int>(bounds[0]);
            // Last bound.
            if (bounds.size() != 1)
              last = convert<int>(bounds[1]);
            else
              last = first;
            string filename;
            // Loop over bounds.
            for (int i = first; i < last + 1; i++)
              {
                // Field base name is empty.
                filename = find_replace(file, "&f", "");
                // Field number is replaced.
                filename = find_replace(filename, "&n", to_str(i));
                // The field name is 'to_str(i)'.
                files_[to_str(i)] = filename;
              }
          }
        else if (vsplit.size() == 1)
          // No numbers.
          {
            // Field base name is replaced.
            file = find_replace(file, "&f", field);
            // No number.
            file = find_replace(file, "&n", "");
            files_[field] = file;
          }
        else
          // With numbers (bounds) and a field base name.
          {
            bounds = split(vsplit[1], "-");
            // First bound.
            first = convert<int>(bounds[0]);
            // Last bound.
            if (bounds.size() != 1)
              last = convert<int>(bounds[1]);
            else
              last = first;
            string filename;
            // Loop over bounds.
            for (int i = first; i < last + 1; i++)
              {
                // Field base name is replaced.
                filename = find_replace(file, "&f", split(field, "_")[0]);
                // Field number is replaced.
                filename = find_replace(filename, "&n", to_str(i));
                // Complete field name is 'vsplit[0] + to_str(i)'.
                files_[vsplit[0] + to_str(i)] = filename;
              }
          }
      } // Loop over all fields and associated files.
  }


  // OTHER METHODS //
  ///////////////////


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_INPUTFILES_CXX
#endif
