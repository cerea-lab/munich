// Copyright (C) 2005-2010, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim, Vivien Mallet
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


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_PHOTOCHEMISTRY_HXX


#include <vector>

#include "AtmoDataHeader.hxx"
#include "BaseModuleParallel.hxx"

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////
  // FORTRAN FUNCTION //
  //////////////////////


#define _chem chem_
#define _chemcl chemcl_
#define _dimensions_racm dimensions_racm_
#define _dimensions_racm2 dimensions_racm2_
#define _dimensions_cb05 dimensions_cb05_
#define _dimensions_cb05_line_ping dimensions_cb05_line_ping_
#define _dimensions_leighton dimensions_leighton_
#define _dimensions_melchior2 dimensions_melchior2_

  extern "C"
  {
    void _chem(int*, int*, int*, int*, int*, int*, double*, double*, double*,
               double*, double*, double*, double*, double*, double*, double*,
               double*, double*, double*, double*, double*, double*, int*,
               double *, double*, double*, int*, double*, double*, int*, int*,
               double*);
    void _chemcl(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
                 double*, double*, double*, double*, double*, double*,
                 double*, double*, double*, double*, double*, double*,
                 double*, double*, double*, int*, double *, double*, double*,
                 double*, double*, double*);

    void _dimensions_racm(int*, int*, int*);
    void _dimensions_racm2(int*, int*, int*);
    void _dimensions_cb05(int*, int*, int*);
    void _dimensions_leighton(int*, int*, int*);
    void _dimensions_melchior2(int*, int*, int*);

  }

  enum CHEMISTRY_MECHANISM
    {
      ZERO,
      RACM,
      RACM2,
      CB05,
      Leighton,
      MELCHIOR2
    };


  ////////////////////
  // PHOTOCHEMISTRY //
  ////////////////////


  /*! \brief This class is a numerical solver for the chemical mechanisms:
    RACM, RACM2, CB05, CB05_line_PinG, Leighton, and MELCHIOR2.
  */
  /*! It uses a second-order Rosenbrock method.
   */
  template<class T>
  class Photochemistry: public BaseModuleParallel
  {

  protected:

    //! Number of species.
    int Ns;
    //! Number of reactions.
    int Nr;
    //! Number of photolysis reactions.
    int Nr_photolysis;
    //! Number of species with volume sources.
    int Ns_source;
    //! Number of levels with volume sources.
    int Nz_source;

    //! Number of sub-cycles.
    int Ncycle;

    //! Sorted list of species names.
    vector<string> species_list;

    //! Molecular weights of species.
    Array<T, 1> molecular_weight;
    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 1> ConversionFactor;
    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 2> ConversionFactorJacobian;

    /*! \brief Map between photolysis reactions names and their indices in
      photolysis reactions.
    */
    map<string, int> photolysis_reaction_name;
    //! Indices of photolysis reactions among other reactions.
    Array<int, 1> photolysis_reaction_index;
    //! Indices of species with volume sources.
    Array<int, 1> source_index;

    //! With adaptive time stepping for gaseous chemistry?
    bool with_adaptive_time_step;
    //! With adaptive time stepping for gaseous chemistry?
    int option_adaptive_time_step;
    T adaptive_time_step_tolerance;
    //! Minimum time step that can be used.
    T min_adaptive_time_step;
    //! Maximum time step that can be used.
    T max_adaptive_time_step;
    /*! \brief Photolysis rates obtained by the tabulation generated by SPACK
      or by binary files from FastJ.
    */
    int option_photolysis_tabulation;
    //! Chemistry mechanism used: RACM, RACM2, CB05, CB05_line_PinG, 
    //!                           Leighton, MELCHIOR2.
    int option_chemistry;


  public:


    /*** Constructor ***/

    Photochemistry();

    /*** Other methods ***/

    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void Forward(ClassModel& Model);

    template<class ClassModel>
    void Backward(ClassModel& Model);

    void Forward(T current_time,
                 Data<T, 3>& Attenuation_i,
                 Data<T, 3>& SpecificHumidity_i,
                 Data<T, 3>& Temperature_i,
                 Data<T, 3>& Pressure_i,
                 Data<T, 4>& Source_i,
                 Data<T, 4>& PhotolysisRate_i,
                 T next_time,
                 Data<T, 3>& Attenuation_f,
                 Data<T, 3>& SpecificHumidity_f,
                 Data<T, 3>& Temperature_f,
                 Data<T, 3>& Pressure_f,
                 Data<T, 4>& Source_f,
                 Data<T, 4>& PhotolysisRate_f,
                 Array<T, 1>& Longitude,
                 Array<T, 1>& Latitude,
                 Data<T, 4>& Concentration);

    void Forward(T current_time,
                 T& attenuation,
                 T& humid,
                 T& temperature,
                 T& pressure,
                 Array<T, 1>& source,
                 Array<T, 1>& photolysis,
                 T next_time,
                 T& attenuation_f,
                 T& humid_f,
                 T& temperature_f,
                 T& pressure_f,
                 Array<T, 1>& source_f,
                 Array<T, 1>& photolysis_f,
                 T& lon, T& lat,
                 Array<T, 1>& concentration);

    void Backward(T current_time,
                  Data<T, 3>& Attenuation_i,
                  Data<T, 3>& SpecificHumidity_i,
                  Data<T, 3>& Temperature_i,
                  Data<T, 3>& Pressure_i,
                  Data<T, 4>& Source_i,
                  Data<T, 4>& PhotolysisRate_i,
                  T next_time,
                  Data<T, 3>& Attenuation_f,
                  Data<T, 3>& SpecificHumidity_f,
                  Data<T, 3>& Temperature_f,
                  Data<T, 3>& Pressure_f,
                  Data<T, 4>& Source_f,
                  Data<T, 4>& PhotolysisRate_f,
                  Array<T, 1>& Longitude,
                  Array<T, 1>& Latitude,
                  Data<T, 4>& Concentration,
                  Data<T, 4>& Source_i_ccl,
                  Data<T, 4>& Source_f_ccl,
                  Data<T, 4>& Concentration_ccl);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_PHOTOCHEMISTRY_HXX
#endif
