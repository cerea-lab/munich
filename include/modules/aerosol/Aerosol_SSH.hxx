// Copyright (C) 2020, CEREA (ENPC - EDF R&D)
// Author(s): Youngseob Kim
//
// This file is a component of the air quality modeling system Polyphemus.
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


#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SSH_HXX

extern "C"
{
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/ipc.h>
}

#include <vector>
#include <map>
#include "AtmoData.hxx"
#include "BaseModuleParallel.cxx"
#include <cmath>
#include <string>
#include <iostream>

#include "API.cxx"
#include <dlfcn.h>

// Global variables.

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////
  // FORTRAN FUNCTION //
  //////////////////////

#ifdef POLYPHEMUS_SINGLE_UNDERSCORE
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#elif defined(__GNUG__) && __GNUG__ < 4 && !defined(__INTEL_COMPILER)
#undef POLYPHEMUS_DOUBLE_UNDERSCORE
#define POLYPHEMUS_DOUBLE_UNDERSCORE
#endif

#define _vsrmchem vsrmchem_
#define _simple_aqueous_module simple_aqueous_module_
  
  extern "C"
  {

    void _vsrmchem(int*, int*, int*, int*, int*,
                   double*, double*, double*, double*, double*,
                   double*, double*, double*, double*, double*,
                   double*, double*, double*, double*, double*,
                   double*, double*, double*,
                   int*,
                   double*, double*,
                   int*, int*, int*,
                   int*, int*, int*, int*, int*);

    void _simple_aqueous_module(int*, int*, int*, int*, int*,
                                double*, double*, double*, double*, double*,
                                double*, double*, double*, double*, double*,
                                double*, double*, double*, double*, double*,
                                double*, double*, double*,
                                int*,
                                double*, double*,
                                int*, int*, int*,
                                int*, int*, int*, int*, int*);
  }


  
  /////////////////
  // AEROSOL_SSH //
  /////////////////


  //! \brief This class is a numerical solver for the chemical mechanism
  //! RACM_SSH, CB05_SSH and RACM2_SSH.
  /*! It uses a second-order Rosenbrock method.
   */
  template<class T>
  class Aerosol_SSH: public BaseModuleParallel
  {

  protected:


    API<T> api;

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

    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 1> ConversionFactor;
    //! Conversion factor from \mu.g/m3 to molecules/cm3.
    Array<T, 2> ConversionFactorJacobian;
	
    // /*! \brief Map between photolysis reactions names and their indices in
    //   photolysis reactions.
    // */
    // map<string, int> photolysis_reaction_name;

    // //! Indices of photolysis reactions among other reactions.
    // Array<int, 1> photolysis_reaction_index;

    //! Indices of species with volume sources.
    Array<int, 1> source_index;

    //! Number of inorganic aerosol species.
    int Ns_inorganic_aer;
   
    //! Number of inert aerosol species.
    int Ns_inert_aer;
    
    //! Number of organic aerosol species.
    int Ns_organic_aer;

    //! Number of aerosol species in the different layers.
    //! N_aerosol_layers in SSH-aerosol    
    int Ns_aer;

    //! Number of aerosol species without considering layers.
    int Ns_aer_nolayer;

    //! Number of aerosol layers.
    int N_layer;    
    
    //! Number of aerosol groups.
    int Ngroup_aer;

    //! Number of aerosol composition sections.
    int Ncomposition_aer;

    //! Number of aerosol fractions.
    int Nfraction_aer;    

    //! Number of aerosol size sections.
    int Nsize_section_aer;
    
    //! Number of aerosol bins.
    int Nbin_aer;
	
    //! Bin index corresponding to fixed cutting diameter.
    int cutting_bin;
    /* bins from 1 to cutting_bin included are at equilibrium
       and are dynamic from cutting_bin+1 to Nbin_aer */

    //! Indices of species with heterogeneous reactions.
    Array<int, 1> heterogeneous_reaction_index;
    
    void *_aerosol_so = NULL;

    T conserving_mass_tolerance;

    // bool is_aerosol_layer;

    bool with_gas_chemistry;

    bool with_external_composition;

    int i_hydrophilic;

    //    vector<string> aerosol_spec_name;
    Array<string, 1> aerosol_spec_name;

    // in g/mol 
    Array<T, 1> ssh_mass_density_layers;

    Array<int, 1> index_groups;
    
    Array<int, 1> index_groups_ext;

    Array<T, 3> discretization_composition;

    Array<T, 3> discretization_composition_conv;
    
    Array<int, 1> aerosol_type;

    Array<T, 1> ssh_diam_input;
    
  public:

    /*** Configuration ***/

    map<string, bool> option_process_aer;
   
    //! Numerical solver for dynamic bin condensation (etr, ros2 or ebi).
    string dynamic_condensation_solver;

    //! Cutting diameter between equilibrium and dynamic bins.
    double fixed_cutting_diameter;

    //! Sulfate condensation computation method (equilibrium, dynamic).
    string sulfate_computation;

    //! SOA computation method (equilibrium, dynamic).
    string soa_computation;

    //! Name of the aqueous module.
    string aqueous_module;

    /*! Redistribution method for aqueous chemistry module.
    */
    string redistribution_method;
    int iredist;
    
    //! Bins bounds (in m).
    Array<T, 1> BinBound_aer;

    //! Aerosol density (kg / m^3).
    T FixedDensity_aer;
    
    //! liquid water content threshold for clouds.
    double lwc_cloud_threshold;


    //! Gas species cloud interacting index ().
    int Ns_cloud_interact;
    vector<int> species_index_cloud_interact;
    
    
    /*** Constructor and destructor ***/

    Aerosol_SSH();
    ~Aerosol_SSH();
    
    /*** Other methods ***/


    template<class ClassModel>
    void InitSharedLib(ClassModel& Model);

    void UpdateSharedLib(T current_time,
                         T attenuation,
                         T humidity,
                         T temperature,
                         T pressure,
                         T delta_t,
                         T lon, T lat);
    
    template<class ClassModel>
    void CheckConfiguration(ClassModel& Model);
    
    template<class ClassModel>
    void Init(ClassModel& Model);

    template<class ClassModel>
    void DisplayConfiguration(ClassModel& Model);
    
    template<class ClassModel>
    void Forward(ClassModel& Model);

    template<class ClassModel>
    void Forward_aer(ClassModel& Model);

    template<class ClassModel>
    void InitDistribution(ClassModel& Model);

    void InitDistribution(Data<T, 4>& Concentration,
                          Data<T, 5>& Concentration_aer,
                          Data<T, 4>& NumberConcentration_aer);

    void InitDistribution(Data<T, 2>& Concentration,
                          Data<T, 3>& Concentration_aer,
                          Data<T, 2>& NumberConcentration_aer);
    
    void InitDistribution(Array<T, 1>& Concentration1D,
                          Array<T, 2>& Concentration_aer2D,
                          Array<T, 1>& NumberConcentration_aer1D,
                          Array<T, 2>& InitConcentration_aer,
                          Array<T, 1>& InitNumberConcentration_aer);

    //! Functions to get informations from SSH-aerosol.
    
    Array<T, 1> GetMassDensityLayers();

    int GetNsAer();
    
    int GetNsAerNoLayer();

    int GetNLayer();

    int GetIHydrophilic();

    int GetNbinAer();

    int GetNSizeSectionAer();    

    bool IsExternalComposition();

    Array<int, 1> GetIndexGroups();

    Array<int, 1> GetIndexGroupsExt();
    
    int GetNgroup_aer();

    int GetNcomposition_aer();

    int GetNfraction_aer();

    Array<T, 3> GetCompositionBounds();

    Array<int, 1> GetAerosolType();

    Array<string, 1> GetAerosolSpecName();
    
    
    void Forward(T current_time,
                 Data<T, 3>& Attenuation_i,
                 Data<T, 3>& SpecificHumidity_i,
                 Data<T, 3>& Temperature_i,
                 Data<T, 3>& Pressure_i,
                 Data<T, 4>& VolumeEmission_i,
                 Data<T, 4>& PhotolysisRate_i,
                 T next_time,
                 Data<T, 3>& Attenuation_f,
                 Data<T, 3>& SpecificHumidity_f,
                 Data<T, 3>& Temperature_f,
                 Data<T, 3>& Pressure_f,
                 Data<T, 4>& VolumeEmission_f,
                 Data<T, 4>& PhotolysisRate_f,
                 Array<T, 1>& Longitude,
                 Array<T, 1>& Latitude,
                 Data<T, 4>& Concentration,
                 Data<T, 3>& LiquidWaterContent_i,
                 Data<T, 4>& WetDiameter_aer,
                 Data<T, 5>& Concentration_aer,
                 Data<T, 3>& pH,
                 Data<T, 4>& NumberConcentration_aer);

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
                 Array<T, 1>& concentration,
                 T& liquid_water_content,
                 Array<T, 1>& wet_diameter_aer,
                 Array<T, 2>& concentration_aer,
                 T& ph,
                 Array<T, 1>& number_concentration_aer);

    // Include the cacul of number concentration for aerosol
    void Forward_aer(T current_time,
                     Data<T, 3>& SpecificHumidity_i,
                     Data<T, 3>& Temperature_i,
                     Data<T, 3>& Pressure_i,
                     T next_time,
                     Data<T, 4>& Concentration,
                     Data<T, 3>& LiquidWaterContent_i,
                     Data<T, 2>& Rain_i,
                     Array<T, 1>& VerticalInterface,
                     Data<T, 5>& Concentration_aer,
                     Data<T, 3>& InCloudWetDepositionFlux,
                     Data<T, 4>& InCloudWetDepositionFlux_aer,
                     Data<T, 3>& pH,
                     Data<T, 4>& NumberConcentration_aer,
                     Data<T, 3>& InCloudWetDepositionFluxNumber_aer);
	
    void Forward_aer(T current_time,
                     T& specifichumidity,
                     T& temperature,
                     T& pressure,
                     T delta_t,
                     Array<T, 1>& concentration,
                     T& liquidwatercontent,
                     T& rain,
                     Array<T, 1>& CurrentVerticalInterface,
                     Array<T, 2>& concentration_aer,
                     Array<T, 1>& incloudwetdepositionflux,
                     Array<T, 2>& incloudwetdepositionflux_aer,
                     T& ph,
                     T& lwc_avg, T& heightfog, int& ifog,
                     Array<T, 1>& number_concentration_aer,
                     Array<T, 1>& incloudwetdepositionfluxnumber_aer,
                     Array<T, 1>& wet_diameter_aer);


    bool IsRequired(string field);
    bool IsComputed(string field);
    void FogSettling(Array<T, 1>& LiquidWaterContent,
                     T& lwc_cloud_threshold,
                     Array<T, 1>& VerticalInterface,
                     T& lwc_avg,
                     int& nfoglay,
                     T& heightfog);

    // void WriteSSHOutput();
    
    // void Read_Coagulation_Coefficient(const string &input_file);
    // void ComputeCoagulationCoefficient(int Nmc);

    
    // Array<ptr_to_real_array, 1> repartition_coefficient;
    // Array<ptr_to_integer_array, 1> index1_repartition_coefficient;
    // Array<ptr_to_integer_array, 1> index2_repartition_coefficient;

    // vector<ptr_to_real_array> repartition_coef;
    
  };


  

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_AEROSOL_AEROSOL_SSH_HXX
#endif
