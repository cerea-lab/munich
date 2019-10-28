// Copyright (C) 2005-2010, ENPC - INRIA - EDF R&D
// Author(s): Youngseob Kim, Vivien Mallet, Meryem Ahmed de Biasi
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


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_PHOTOCHEMISTRY_CXX


#include "Photochemistry.hxx"

#include "AtmoData.hxx"
#include "BaseModuleParallel.cxx"


namespace Polyphemus
{


  //! Default constructor.
  template<class T>
  Photochemistry<T>::Photochemistry()
  {
  }


  //! Initialization of the scheme.
  /*! In case the model is incompatible with the chemical mechanism, an
    exception is thrown.
    \param Model model with the following interface:
    <ul>
    <li> GetNs()
    <li> GetSpeciesList()
    <li> GetSpeciesFile()
    <li> GetNr_photolysis()
    <li> GetPhotolysisReactionList()
    <li> GetNs_source()
    <li> GetNz_source()
    <li> SourceGlobalIndex(int)
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Photochemistry<T>::Init(ClassModel& Model)
  {
    Ncycle = 1;

    ConfigStream config_main(Model.GetConfigurationFile());

    config_main.SetSection("[options]");

    string str_option_chemistry;
    config_main.PeekValue("Option_chemistry", str_option_chemistry);

    config_main.PeekValue("With_adaptive_time_step_for_gas_chemistry",
                          with_adaptive_time_step);
    config_main.PeekValue("Photolysis_tabulation_option", option_photolysis_tabulation);

    // Adaptive time step.
    if (with_adaptive_time_step)
      {
        option_adaptive_time_step = 1;
        config_main.PeekValue("Adaptive_time_step_tolerance",
                              adaptive_time_step_tolerance);
        config_main.PeekValue("Min_adaptive_time_step",
                              min_adaptive_time_step);
        config_main.PeekValue("Max_adaptive_time_step",
                              max_adaptive_time_step);
      }
    else
      option_adaptive_time_step = 0;

    // Option for chemical mechanisms.
    if (str_option_chemistry == "RACM")
      {
        option_chemistry = RACM;
        _dimensions_racm(&Ns, &Nr, &Nr_photolysis);
      }
    else if (str_option_chemistry == "RACM2")
      {
        option_chemistry = RACM2;
        _dimensions_racm2(&Ns, &Nr, &Nr_photolysis);
      }
    else if (str_option_chemistry == "CB05")
      {
        option_chemistry = CB05;
        _dimensions_cb05(&Ns, &Nr, &Nr_photolysis);
      }
    else if (str_option_chemistry == "Leighton")
      {
     	option_chemistry = Leighton;
	    _dimensions_leighton(&Ns,&Nr,&Nr_photolysis);
      }
    else
      throw string("Wrong index for option_chemistry \n");

    if (Model.GetNs() != Ns)
      throw string("Incompatibility: model manages ") + to_str(Model.GetNs())
        + string(" species while chemistry has ") + to_str(Ns) + " species.";
    if (Model.GetNr_photolysis() != Nr_photolysis)
      throw string("Incompatibility: model manages ")
        + to_str(Model.GetNr_photolysis())
        + string(" photolysis reactions while chemistry has ")
        + to_str(Nr_photolysis) + " photolysis reactions.";
    species_list = Model.GetSpeciesList();

    photolysis_reaction_index.resize(Nr_photolysis);
    molecular_weight.resize(Ns);

    int i, j;

    const string& config_file = Model.GetConfigurationFile();
    const string& species_file = Model.GetSpeciesFile();
    ConfigStreams config(config_file, species_file);

    string molecular_weight_file = species_file;
    config.SetSection("[molecular_weight]");
    // Checks if there is an included file.
    if (config.GetElement() == "file")
      molecular_weight_file = config.GetElement();
    ConfigStream molecular_weight_config(molecular_weight_file);
    molecular_weight_config.SetSection("[molecular_weight]");
    for (int i = 0; i < Ns; i++)
      molecular_weight_config.PeekValue(species_list[i], molecular_weight(i));

    config.SetSection("[photolysis_reaction_index]");
    string species;
    for (i = 0; i < Nr_photolysis; i++)
      {
        species = config.GetElement();
        photolysis_reaction_name[species] = convert<int>(config.GetElement());
      }

    for (i = 0; i < Nr_photolysis; i++)
      photolysis_reaction_index(i) =
        photolysis_reaction_name[Model.GetPhotolysisReactionList()[i]];
    Ns_source = 0;

    // Conversions.
    double Navogadro = 6.02213e23;
    ConversionFactor.resize(Ns);
    for (i = 0; i < Ns; i++)
      ConversionFactor(i) = Navogadro * 1e-12 / molecular_weight(i);
    ConversionFactorJacobian.resize(Ns, Ns);
    for (i = 0; i < Ns; i++)
      for (j = 0; j < Ns; j++)
        ConversionFactorJacobian(i, j) = molecular_weight(j)
          / molecular_weight(i);

    Ns_source = Model.GetNs_source();
    Nz_source = Model.GetNz_source();
    source_index.resize(Ns_source);
    for (int i = 0; i < Ns_source; i++)
      source_index(i) = Model.SourceGlobalIndex(i);

    // To be called even if there is no parallelization.
    BaseModuleParallel::Init(Model);

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Partition along axis x.
    BuildPartition_x();
#endif
  }


  //! Performs an integration over one time step.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("Attenuation_i
    <li> D3("SpecificHumidity_i"),
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetSource_i()
    <li> D4("PhotolysisRate_i")
    <li> D3("Attenuation_f")
    <li> D3("SpecificHumidity_f")
    <li> D3("Temperature_f")
    <li> D3("Pressure_f")
    <li> GetSource_f()
    <li> D4("PhotolysisRate_f")
    <li> GetGridXArray1D()
    <li> GetGridYArray1D()
    <li> GetConcentration()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Photochemistry<T>::Forward(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    if (Model.source_splitting)
      {
        ScatterSlice_x_MPI(Model.GetSource_i());
        ScatterSlice_x_MPI(Model.GetSource_f());
      }
#endif

    Forward(T(date_i.GetNumberOfSeconds()), Model.D3("Attenuation_i"),
            Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
            Model.D3("Pressure_i"), Model.GetSource_i(),
            Model.D4("PhotolysisRate_i"), T(date_f.GetSecondsFrom(date_i)),
            Model.D3("Attenuation_f"), Model.D3("SpecificHumidity_f"),
            Model.D3("Temperature_f"), Model.D3("Pressure_f"),
            Model.GetSource_f(), Model.D4("PhotolysisRate_f"),
            Model.GetGridXArray1D(), Model.GetGridYArray1D(),
            Model.GetConcentration());
  }


  //! Performs one step of backward integration of adjoint model.
  /*!
    \param Model model with the following interface:
    <ul>
    <li> GetCurrentDate()
    <li> GetNextDate()
    <li> D3("Attenuation_i
    <li> D3("SpecificHumidity_i"),
    <li> D3("Temperature_i")
    <li> D3("Pressure_i")
    <li> GetSource_i()
    <li> D4("PhotolysisRate_i")
    <li> D3("Attenuation_f")
    <li> D3("SpecificHumidity_f")
    <li> D3("Temperature_f")
    <li> D3("Pressure_f")
    <li> GetSource_f()
    <li> D4("PhotolysisRate_f")
    <li> GetGridXArray1D()
    <li> GetGridYArray1D()
    <li> GetConcentration()
    <li> D4("Source_i_ccl")
    <li> D4("Source_f_ccl")
    <li> GetConcentration_ccl()
    </ul>
  */
  template<class T>
  template<class ClassModel>
  void Photochemistry<T>::Backward(ClassModel& Model)
  {
    Date date_i = Model.GetCurrentDate();
    Date date_f = Model.GetNextDate();

    Backward(T(date_i.GetNumberOfSeconds()), Model.D3("Attenuation_i"),
             Model.D3("SpecificHumidity_i"), Model.D3("Temperature_i"),
             Model.D3("Pressure_i"), Model.GetSource_i(),
             Model.D4("PhotolysisRate_i"), T(date_f.GetSecondsFrom(date_i)),
             Model.D3("Attenuation_f"), Model.D3("SpecificHumidity_f"),
             Model.D3("Temperature_f"), Model.D3("Pressure_f"),
             Model.GetSource_f(), Model.D4("PhotolysisRate_f"),
             Model.GetGridXArray1D(), Model.GetGridYArray1D(),
             Model.GetConcentration(), Model.D4("Source_i_ccl"),
             Model.D4("Source_f_ccl"), Model.GetConcentration_ccl());
  }


  //! Performs an integration over one time step.
  /*!
    \param current_time starting time in seconds.
    \param Attenuation_i cloud attenuation coefficients at the beginning of
    the time step.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param Source_i volume sources at the beginning of the time step.
    \param PhotolysisRate_i photolysis rates at the beginning of the time
    step.
    \param Attenuation_f cloud attenuation coefficients at the end of the
    time step.
    \param SpecificHumidity_f specific humidity at the end of the time step.
    \param Temperature_f temperature at the end of the time step.
    \param Pressure_f pressure at the end of the time step.
    \param Source_f volume sources at the end of the time step.
    \param PhotolysisRate_f photolysis rates at the end of the time step.
    \param Longitude longitudes.
    \param Latitude latitudes.
    \param Concentration concentrations.
  */
  template<class T>
  void Photochemistry<T>::Forward(T current_time,
                                  Data<T, 3>& Attenuation_i,
                                  Data<T, 3>& SpecificHumidity_i,
                                  Data<T, 3>& Temperature_i,
                                  Data<T, 3>& Pressure_i,
                                  Data<T, 4>& Source_i,
                                  Data<T, 4>& PhotolysisRate_i,
                                  T delta_t,
                                  Data<T, 3>& Attenuation_f,
                                  Data<T, 3>& SpecificHumidity_f,
                                  Data<T, 3>& Temperature_f,
                                  Data<T, 3>& Pressure_f,
                                  Data<T, 4>& Source_f,
                                  Data<T, 4>& PhotolysisRate_f,
                                  Array<T, 1>& Longitude,
                                  Array<T, 1>& Latitude,
                                  Data<T, 4>& Concentration)
  {
    int first_index_along_x, last_index_along_x;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // This array will store the concentration data with the species
    // axis permuted with the X-axis.
    Array<T, 4> OrdConc;
    ScatterSlice_x_MPI(Concentration, OrdConc);
    GetEdgePartition_x(first_index_along_x, last_index_along_x);
#else
    first_index_along_x = 0;
    last_index_along_x = Concentration.GetLength(3);
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp = GetNthreads_openmp();
#pragma omp parallel for num_threads(Nthreads_openmp)	\
  firstprivate(first_index_along_x, last_index_along_x)
#endif
    for (int i = first_index_along_x; i < last_index_along_x; i++)
      for (int k = 0; k < Concentration.GetLength(1); k++)
        for (int j = 0; j < Concentration.GetLength(2); j++)
          {
            Array<T, 1> Photolysis1D(Nr_photolysis);
            Array<T, 1> Photolysis1D_f(Nr_photolysis);
            Array<T, 1> Concentration1D(Ns);
            Array<T, 1> Source1D(Ns_source);
            Array<T, 1> Source1D_f(Ns_source);
            int s;

            for (s = 0; s < Nr_photolysis; s++)
              {
                Photolysis1D(s) = PhotolysisRate_i()(s, k, j, i);
                Photolysis1D_f(s) = PhotolysisRate_f()(s, k, j, i);
              }

            for (s = 0; s < Ns; s++)
              {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
                // Rank of indices must match the one defined in
                // ScatterSlice_x_MPI.
                Concentration1D(s) = OrdConc(i, k, j, s);
#else
                Concentration1D(s) = Concentration(s, k, j, i);
#endif
              }

            if (k < Nz_source)
              {
                for (s = 0; s < Ns_source ; s++)
                  {
                    Source1D(s) = Source_i()(s, k, j, i);
                    Source1D_f(s) = Source_f()(s, k, j, i);
                  }
              }
            else
              {
                Source1D = 0;
                Source1D_f = 0;
              }

            Forward(current_time, Attenuation_i(k, j, i),
                    SpecificHumidity_i(k, j, i), Temperature_i(k, j, i),
                    Pressure_i(k, j, i), Source1D, Photolysis1D,
                    delta_t, Attenuation_f(k, j, i),
                    SpecificHumidity_f(k, j, i), Temperature_f(k, j, i),
                    Pressure_f(k, j, i), Source1D_f, Photolysis1D_f,
                    Longitude(i), Latitude(j), Concentration1D);

            for (s = 0; s < Ns; s++)
              {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
                // Rank of indices must match the one defined in
                // ScatterSlice_x_MPI.
                OrdConc(i, k, j, s) = Concentration1D(s);
#else
                Concentration(s, k, j, i) = Concentration1D(s);
#endif
              }
          }
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    GatherSlice_x_MPI(OrdConc, Concentration);
#endif
  }


  //! Performs an integration over one time step at one location.
  /*!
    \param current_time starting time in seconds.
    \param Attenuation_i cloud attenuation coefficients at the beginning of
    the time step.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param Source_i volume sources at the beginning of the time step.
    \param PhotolysisRate_i photolysis rates at the beginning of the time
    step.
    \param Attenuation_f cloud attenuation coefficients at the end of the
    time step.
    \param SpecificHumidity_f specific humidity at the end of the time step.
    \param Temperature_f temperature at the end of the time step.
    \param Pressure_f pressure at the end of the time step.
    \param Source_f volume sources at the end of the time step.
    \param PhotolysisRate_f photolysis rates at the end of the time step.
    \param Longitude longitudes.
    \param Latitude latitudes.
    \param Concentration concentrations.
  */
  template<class T>
  void Photochemistry<T>::Forward(T current_time,
                                  T& attenuation,
                                  T& humid,
                                  T& temperature,
                                  T& pressure,
                                  Array<T, 1>& source,
                                  Array<T, 1>& photolysis,
                                  T delta_t,
                                  T& attenuation_f,
                                  T& humid_f,
                                  T& temperature_f,
                                  T& pressure_f,
                                  Array<T, 1>& source_f,
                                  Array<T, 1>& photolysis_f,
                                  T& lon, T& lat,
                                  Array<T, 1>& concentration)
  {
    _chem(&Ns, &Nr, &Nr_photolysis, photolysis_reaction_index.data(),
          &Ns_source, source_index.data(),
          ConversionFactor.data(), ConversionFactorJacobian.data(),
          &current_time, &attenuation,
          &humid, &temperature,
          &pressure, source.data(),
          photolysis.data(), &delta_t, &attenuation_f,
          &humid_f, &temperature_f,
          &pressure_f, source_f.data(),
          photolysis_f.data(), &Ncycle, &lon, &lat,
          concentration.data(),
          &option_adaptive_time_step, &adaptive_time_step_tolerance,
          &min_adaptive_time_step,
          &option_photolysis_tabulation, &option_chemistry, &max_adaptive_time_step);
  }


  //! Performs one step of backward integration of adjoint model.
  /*!
    \param current_time starting time in seconds.
    \param Attenuation_i cloud attenuation coefficients at the beginning of
    the time step.
    \param SpecificHumidity_i specific humidity at the beginning of the time
    step.
    \param Temperature_i temperature at the beginning of the time step.
    \param Pressure_i pressure at the beginning of the time step.
    \param Source_i volume sources at the beginning of the time step.
    \param PhotolysisRate_i photolysis rates at the beginning of the time
    step.
    \param Attenuation_f cloud attenuation coefficients at the end of the
    time step.
    \param SpecificHumidity_f specific humidity at the end of the time step.
    \param Temperature_f temperature at the end of the time step.
    \param Pressure_f pressure at the end of the time step.
    \param Source_f volume sources at the end of the time step.
    \param PhotolysisRate_f photolysis rates at the end of the time step.
    \param Longitude longitudes.
    \param Latitude latitudes.
    \param Concentration concentrations.
    \param Source_i_ccl adjoint volume sources data at the beginning of the
    time step.
    \param Source_f_ccl adjoint volume sources data at the end of the
    time step.
    \param Concentration_ccl adjoint concentration data.
  */
  template<class T>
  void Photochemistry<T>::Backward(T current_time,
                                   Data<T, 3>& Attenuation_i,
                                   Data<T, 3>& SpecificHumidity_i,
                                   Data<T, 3>& Temperature_i,
                                   Data<T, 3>& Pressure_i,
                                   Data<T, 4>& Source_i,
                                   Data<T, 4>& PhotolysisRate_i,
                                   T delta_t,
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
                                   Data<T, 4>& Concentration_ccl)
  {
    int Nz = Concentration.GetLength(1);
    int Ny = Concentration.GetLength(2);
    int Nx = Concentration.GetLength(3);

    _chemcl(&Nx, &Ny, &Nz, &Ns, &Nr,
            &Nr_photolysis, photolysis_reaction_index.data(),
            &Ns_source, source_index.data(), &Nz_source,
            molecular_weight.data(),
            &current_time, Attenuation_i.GetData(),
            SpecificHumidity_i.GetData(), Temperature_i.GetData(),
            Pressure_i.GetData(), Source_i.GetData(),
            PhotolysisRate_i.GetData(), &delta_t, Attenuation_f.GetData(),
            SpecificHumidity_f.GetData(), Temperature_f.GetData(),
            Pressure_f.GetData(), Source_f.GetData(),
            PhotolysisRate_f.GetData(), &Ncycle, Longitude.data(),
            Latitude.data(), Concentration.GetData(),
            Source_i_ccl.GetData(), Source_f_ccl.GetData(),
            Concentration_ccl.GetData());
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_PHOTOCHEMISTRY_CXX
#endif
