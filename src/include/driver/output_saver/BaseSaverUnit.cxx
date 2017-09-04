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


#ifndef POLYPHEMUS_FILE_OUTPUT_SAVER_BASESAVERUNIT_CXX


#include "BaseSaverUnit.hxx"


//! Checks equality with a tolerance interval of epsilon.
/*!
  \param x first number.
  \param y second number.
  \param epsilon relative tolerance.
  \return True if \a x equals \a y with a relative tolerance of \a epsilon.
*/
template<class T>
bool is_equal(T x, T y, T epsilon = 1.e-6)
{
  return abs(x - y) <=  0.5 * epsilon * (abs(x) + abs(y));
}


//! Checks whether a number is multiple of another, with given tolerance.
/*!
  \param x possible multiple.
  \param d base number.
  \param epsilon relative tolerance.
  \return True if \a x is a multiple of \a d with a relative tolerance of
  epsilon.
*/
template<class T>
bool is_multiple(T x, T d, T epsilon = 1.e-6)
{
  int i = int(x / d + .5);
  return is_equal(x, T(i) * d, epsilon);
}


namespace Polyphemus
{


  ///////////////////
  // BASESAVERUNIT //
  ///////////////////


  //! Main constructor.
  template<class T, class ClassModel>
  BaseSaverUnit<T, ClassModel>::BaseSaverUnit()
  {
  }


  //! Destructor.
  template<class T, class ClassModel>
  BaseSaverUnit<T, ClassModel>::~BaseSaverUnit()
  {
  }


  //! Returns the group of the saver unit.
  /*!
    \return The group of the saver unit, that is, "forecast".
  */
  template<class T, class ClassModel>
  string BaseSaverUnit<T, ClassModel>::GetGroup() const
  {
    return "forecast";
  }


  //! First initialization.
  /*! Reads the configuration.
    \param config_stream configuration stream.
    \param Model model with the following interface:
    <ul>
    <li> GetX_min()
    <li> GetDelta_x()
    <li> GetNx()
    <li> GetY_min()
    <li> GetDelta_y()
    <li> GetNy()
    <li> GetNz()
    <li> GetCurrentDate()
    <li> GetDelta_t()
    <li> GetNt()
    </ul>
  */
  template<class T, class ClassModel>
  void BaseSaverUnit<T, ClassModel>::Init(ConfigStream& config,
					  ClassModel& Model)
  {
    // Initial and final dates for saving.
    Date Date_min = Model.GetCurrentDate();
    Date Date_max(Date_min);
    Date_max.AddSeconds(Model.GetNt() * Model.GetDelta_t());

    string element;
    config.PeekValue("Date_beg", element);
    if (is_num(element) && to_num<int>(element) == -1)
      date_beg = Date_min;
    else
      date_beg = config.PeekValue("Date_beg");

    config.PeekValue("Date_end", element);
    if (is_num(element) && to_num<int>(element) == -1)
      date_end = Date_max;
    else
      date_end = config.PeekValue("Date_end");

    if (date_beg > date_end || date_beg < Date_min || date_end > Date_max)
      throw string("Saving dates are not inside the simulation period.");
      
    if (!is_multiple(T(date_beg.GetSecondsFrom(Date_min)),
                     Model.GetDelta_t()))
      throw string("\"Date_beg\" does not match a simulation timestep.");

    if (!is_multiple(T(Date_max.GetSecondsFrom(date_end)),
                     Model.GetDelta_t()))
      throw string("\"Date_end\" does not match a simulation timestep.");
      
    // Interval length.
    config.PeekValue("Interval_length", element);
    if (is_num(element))
      interval_length = to_num<int>(element);
    else
      {
	T Delta_t_save;
	element = lower_case(element);
	if (element == "hourly")
	  Delta_t_save = 3600;
	else if (element == "daily")
	  Delta_t_save = 86400;
	else
	  throw string("Value of \"Interval_length\" not recognized.");
	
	T Delta_t_model = Model.GetDelta_t();
	if (is_multiple(Delta_t_save, Delta_t_model))
	  interval_length = int(Delta_t_save / Delta_t_model + .5);
	else
	  throw string("Value of \"Interval_length\" must be a")
	    + string(" multiple of model timestep.");
      }

    if (GetType() != "nesting" && GetType() != "nesting_aer")
      {
	config.PeekValue("Averaged", averaged);
	config.PeekValue("Initial_concentration", initial_concentration);
      }

    config.Find("Species");
    species_list = split(config.GetLine());
    if (GetType() != "domain_aer"  && GetType() != "dry_deposition_aer"
	&& GetType() != "wet_deposition_aer" && GetType() != "nesting_aer"
	&& GetType() != "subdomain_aer" && GetType() != "indices_list_aer"
	&& GetType() != "coordinates_list_aer" && species_list[0] == "all")
      {
	species_list.clear();
	if (GetType() == "dry_deposition")
	  species_list = Model.GetSpeciesList("DepositionVelocity");
	else if (GetType() == "wet_deposition")
	  species_list = Model.GetSpeciesList("ScavengingCoefficient");
	else
	  species_list = Model.GetSpeciesList();
      }

    // Dimensions of the underlying model.
    base_x_min = Model.GetX_min();
    base_Delta_x = Model.GetDelta_x();
    base_Nx = Model.GetNx();
    base_y_min = Model.GetY_min();
    base_Delta_y = Model.GetDelta_y();
    base_Ny = Model.GetNy();
    base_Nz = Model.GetNz();

    counter = 0;
  }


  //! Initializes the saver at the beginning of each step.
  /*!
    \param Model model (dummy argument).
  */
  template<class T, class ClassModel>
  void BaseSaverUnit<T, ClassModel>::InitStep(ClassModel& Model)
  {
    counter++;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_OUTPUT_SAVER_BASESAVERUNIT_CXX
#endif
