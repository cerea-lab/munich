// Copyright (C) 2005-2016, ENPC - INRIA - EDF R&D
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


#ifndef POLYPHEMUS_FILE_MODULES_CHEMISTRY_PHOTOCHEMISTRYSING_CXX


#include "PhotochemistrySinG.hxx"


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
    config_main.PeekValue("Photolysis_option",option_photolysis);

    // Adaptive time step.
    if (with_adaptive_time_step)
      {
	option_adaptive_time_step = 1;
	config_main.PeekValue("Adaptive_time_step_tolerance",
                              adaptive_time_step_tolerance);
	config_main.PeekValue("Min_adaptive_time_step",
                              min_adaptive_time_step);
      }
    else
      option_adaptive_time_step = 0;

    // Option for chemical mechanisms.
    if (str_option_chemistry == "RACM")
      {
	option_chemistry = RACM;
	_dimensions_racm(&Ns,&Nr,&Nr_photolysis);
	cout<< "Chemical mechanism: RACM \n";
      }
    else if (str_option_chemistry == "RACM2")
      {
	option_chemistry = RACM2;
	_dimensions_racm2(&Ns,&Nr,&Nr_photolysis);
	cout<< "Chemical mechanism: RACM2 \n";
      }
    else if (str_option_chemistry == "CB05")
      {
	option_chemistry = CB05;
	_dimensions_cb05(&Ns,&Nr,&Nr_photolysis);
	cout<< "Chemical mechanism: CB05 \n";
      }
    else if (str_option_chemistry == "Leighton")
      {
	option_chemistry = Leighton;
	_dimensions_leighton(&Ns,&Nr,&Nr_photolysis);
	cout<< "Chemical mechanism: Leighton \n";
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

    ConfigStream config(Model.GetSpeciesFile());
    config.SetSection("[molecular_weight]");
    int i, j;
    for (i = 0; i < Ns; i++)
      config.PeekValue(Model.GetSpeciesList()[i], molecular_weight(i));

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
	  &option_photolysis, &option_chemistry);
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_CHEMISTRY_PHOTOCHEMISTRYSING_CXX
#endif
