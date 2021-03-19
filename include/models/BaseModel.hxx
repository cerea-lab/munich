// Copyright (C) 2005-2018, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet, Shupeng Zhu
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


#ifndef POLYPHEMUS_FILE_MODELS_BASEMODEL_HXX


#include "SeldonDataHeader.hxx"
//#include "BasePointEmission.hxx"
#include "InputFiles.hxx"
#include <vector>
#include <map>
#include <cmath>


namespace Polyphemus
{


  using namespace std;
  using namespace SeldonData;
  using namespace AtmoData;
  using Talos::to_num;
  using Talos::to_str;


  ///////////////
  // BASEMODEL //
  ///////////////


  /*! \brief 'BaseModel' collects the attributes and the methods that are
    shared by most models. */
  /*! 'BaseModels' mainly collects the grids, several option keys, and methods
    to read data on file and to interpolate them in time.
  */
  template<class T>
  class BaseModel
  {

  public:

    /*** Type declarations ***/

    typedef typename map<string, InputFiles<T> >
    ::iterator input_files_iterator;

  protected:

    /*** Configuration ***/

    //! Name of the main configuration file.
    string file_config;
    //! Configuration stream for the main configuration file.
    ConfigStream config;

    //! Options to enable or disable physical processes.
    map<string, bool> option_process;
    /*! \brief Variables that are managed (that is, read or computed) by the
      model or not. */
    map<string, bool> option_manage;

    //! List of input-files groups.
    map<string, InputFiles<T> > input_files;

    /*** Domain ***/

    //! Simulation starting date.
    Date Date_min;
    //! Time step in seconds.
    T Delta_t;
    //! Number of time steps.
    int Nt;

    //! File storing species list.
    string file_species;
    //! Species list.
    vector<string> species_list;
    //! Number of species.
    int Ns;
    //! List of aerosol species.
    vector<string> species_list_aer;
    //! List of aerosol related variables (species mass + number).
    vector<string> list_aer;
    //! Number of aerosol species.
    int Ns_aer;
    //! List of aerosol groups.
    vector<string> groups_list_aer;
    //! Number of aerosol groups.
    int Ngroup_aer;
    //! Number of bins (for aerosols).
    //!(internal: equal to Nsize_section_aer;
    //!external: equal to Nsize_section_aer*Ncomposition_aer)
    int Nbin_aer;
    //! Number of size sections (for aerosols).
    int Nsize_section_aer;
    //! Number of aerosol compositions
    int Ncomposition_aer;

    //! File storing vertical-levels coordinates.
    string file_vertical_levels;
    //! Number of vertical layers.
    int Nz;
    //! Vertical levels.
    Array<T, 1> LayerInterface;
    //! Domain origin along y.
    T y_min;
    //! Space step along y.
    T Delta_y;
    //! Number of points along y.
    int Ny;
    //! Domain origin along x.
    T x_min;
    //! Space step along x.
    T Delta_x;
    //! Number of points along x.
    int Nx;

    //! 3D grid for species.
    RegularGrid<T> GridS3D;
    //! 4D grid for species.
    RegularGrid<T> GridS4D;

    //! 3D grid along z.
    RegularGrid<T> GridZ3D;
    //! 3D grid for interfaces along z.
    RegularGrid<T> GridZ3D_interf;
    //! 4D grid along z.
    RegularGrid<T> GridZ4D;
    //! 4D grid for interfaces along z.
    RegularGrid<T> GridZ4D_interf;

    //! 2D grid along y.
    RegularGrid<T> GridY2D;
    //! 3D grid along y.
    RegularGrid<T> GridY3D;
    //! 3D grid for interfaces along y.
    RegularGrid<T> GridY3D_interf;
    //! 4D grid along y.
    RegularGrid<T> GridY4D;
    //! 4D grid for interfaces along y.
    RegularGrid<T> GridY4D_interf;

    //! 2D grid along x.
    RegularGrid<T> GridX2D;
    //! 3D grid along x.
    RegularGrid<T> GridX3D;
    //! 3D grid for interfaces along x.
    RegularGrid<T> GridX3D_interf;
    //! 4D grid along x.
    RegularGrid<T> GridX4D;
    //! 4D grid for interfaces along x.
    RegularGrid<T> GridX4D_interf;

    /*** Integration ***/

    //! Date of the previous time-step.
    Date previous_date;
    //! Date of the current time-step.
    Date current_date;
    //! Date of the next time-step.
    Date next_date;
    //! Date at which the data currently is.
    Date data_date;
    //! Has the data already been processed in InitStep?
    bool data_gone_through_initstep;
    /*! \brief Has the data already been transformed in Forward? In Forward,
      input data can sometimes be transformed, after which this flag is
      activated. */
    bool data_gone_through_forward;
    //! Number of seconds between the simulation beginning and the current
    //! date.
    T current_time;
    //! Number of steps achieved so far.
    int step;
    //! Is backward integration? Initially false.
    bool backward;

    /*** Fields ***/

    //! Map from 2D field names to their arrays.
    map<string, Array<T, 2>* > A2_map;
    //! Map from 2D field names to their data.
    map<string, Data<T, 2>* > D2_map;
    //! Map from 3D field names to their arrays.
    map<string, Array<T, 3>* > A3_map;
    //! Map from 3D field names to their data.
    map<string, Data<T, 3>* > D3_map;
    //! Map from 4D field names to their arrays.
    map<string, Array<T, 4>* > A4_map;
    //! Map from 4D field names to their data.
    map<string, Data<T, 4>* > D4_map;
    //! Map from 5D field names to their arrays.
    map<string, Array<T, 5>* > A5_map;
    //! Map from 5D field names to their data.
    map<string, Data<T, 5>* > D5_map;

    //! Map from field names to species lists.
    map<string, vector<string>* > field_species;

    //! Map from field names to bins lists.
    map<string, vector<int>* > field_bins;

    /*! List of input parameters that are considered input of the model time
      stepper, except the state vector. */
    vector<string> parameter_name;

    // /*** Source terms ***/

    // //! Interface to manage point emissions.
    // BasePointEmission<T>* PointEmissionManager; // YK

    /*** State ***/

    //! Concentrations.
    Data<T, 4> Concentration;
    //! Adjoint data for concentrations.
    Data<T, 4> Concentration_ccl;
    Data<T, 3> Vertical_Wind;
    //! Aerosols concentrations.
    Data<T, 5> Concentration_aer;
	//! Aerosols number concentrations.
    Data<T, 4> NumberConcentration_aer;

    // //**LL: Remove stationary regime
    // //! Time step chemistry in seconds.
    // T Delta_t_chem;
    //LL*
    // string file_species_aer;
    // vector<string> groups_list_aer;
    // vector<string> list_aer;
    // int Nsize_section_aer;
    // int Ngroup_aer;
    // int Ncomposition_aer;
    string scavenging_model;
    T building_density;
    //LL
    
  public:

    /*** Constructors and destructor ***/

    BaseModel();
    BaseModel(string config_file);
    void Construct(string config_file);
    virtual ~BaseModel();

    /*** Configuration ***/

    virtual void ReadConfiguration();
    virtual void CheckConfiguration();

    /*** Initializations ***/

    virtual void Allocate();
    virtual void Init();
    virtual void InitStep();

    virtual void SetDate(Date date);
    virtual void StepBack();
    virtual void StepBack(const Array<T, 4>& concentration);

    /*** Integration ***/

    virtual void Forward();
    bool HasFinished() const;
    virtual void SetBackward(bool flag);
    virtual void Backward();

    /*** Access methods ***/

    string GetConfigurationFile() const;

    bool IsBackward();
    bool HasDataGoneThroughInitStep();
    bool HasDataGoneThroughInitStep(bool new_value);
    bool HasDataGoneThroughForward();
    bool HasDataGoneThroughForward(bool new_value);

    T GetDelta_t() const;
    int GetNt() const;

    string GetSpeciesFile() const;
    int GetNs() const;
    int GetNumSpecies() const;
    int GetNs_aer() const;
    int GetNumSpecies_aer() const;
    int GetNbin_aer() const;
    int GetNgroup_aer() const;
    int GetNcomposition_aer() const;
    int GetNsize_section_aer() const;
    
    Date GetDate_min() const;
    int GetNz() const;
    T GetY_min() const;
    T GetDelta_y() const;
    int GetNy() const;
    T GetX_min() const;
    T GetDelta_x() const;
    int GetNx() const;

    vector<string> GetSpeciesList() const;
    string GetSpeciesName(int i) const;
    bool IsSpecies(string name) const;
    int GetSpeciesIndex(string species) const;
    int GetSpeciesIndex(string species,
                        const vector<string>& ref_species_list) const;
    int GetSpeciesIndex(string field, string species);
    
    vector<string> GetGroupList_aer() const;
    vector<string> GetSpeciesList_aer() const;
    string GetGroupName_aer(int i) const;
    string GetSpeciesName_aer(int i) const;
	string GetName_aer(int i) const;
    bool IsSpecies_aer(string name) const;
    int GetSpeciesIndex_aer(string species) const;
    int GetGroupIndex_aer(string species) const;
        
    const vector<string>& GetSpeciesList(string field);
    const vector<int>& GetBinsList_aer(string field);

    Array<T, 1>& GetGridZArray1D();
    Array<T, 1>& GetGridYArray1D();
    Array<T, 1>& GetGridXArray1D();
    Array<T, 1>& GetLayerInterface();

    Date GetPreviousDate() const;
    Date GetCurrentDate() const;
    T GetCurrentTime() const;
    Date GetNextDate() const;

    bool HasField(string field);
    string FindField(string field);
    int GetFieldDimension(string field);
    Array<T, 2>& A2(string field);
    Data<T, 2>& D2(string field);
    Array<T, 3>& A3(string field);
    Data<T, 3>& D3(string field);
    Array<T, 4>& A4(string field);
    Data<T, 4>& D4(string field);
    Array<T, 5>& A5(string field);
    Data<T, 5>& D5(string field);

    void RegisterParameter(string name);
    int GetNparameter() const;
    string GetParameterName(int i) const;

    // BasePointEmission<T>* GetPointEmission(); // YK

    const Data<T, 4>& GetConcentration() const;
    Data<T, 4>& GetConcentration();
    const Data<T, 4>& GetConcentration_ccl() const;
    Data<T, 4>& GetConcentration_ccl();
    const Data<T, 5>& GetConcentration_aer() const;
    Data<T, 5>& GetConcentration_aer();
    Data<T, 3>& GetVertical_Wind();
    const Data<T, 4>& GetNumberConcentration_aer() const;
    Data<T, 4>& GetNumberConcentration_aer();
    virtual T GetConcentration(int species, T z, T y, T x);
    virtual T GetIntegratedConcentration(int species, T z, T y, T x,
                                         T lz, T ly, T lx);
    virtual T GetConcentration_aer(int species, int diameter, T z, T y, T x);
    virtual T GetNumberConcentration_aer(int diameter, T z, T y, T x);
    virtual bool HasNumberConcentration_aer();

    virtual void ComputeConcentration(const vector<int>& species_index,
                                      const vector<int>& levels);
    virtual void ComputeConcentration();

    bool IsManaged(string field) const;
    void Manage(string field, bool status = true);

    /*** Time managements ***/

    void SetCurrentDate(Date date);
    void AddTime(T time);
    void SubtractTime(T time);

  protected:

    void Register(Data<T, 2>& Data_i, Data<T, 2>& Data_f, string field_name);
    void Register(Data<T, 2>& Data_i, string field_name);

    void Register(Data<T, 3>& Data_i, Data<T, 3>& Data_f, string field_name);
    void Register(Data<T, 3>& Data_i, string field_name);

    void Register(Data<T, 4>& Data_i, Data<T, 4>& Data_f, string field_name);
    void Register(Data<T, 4>& Data_i, string field_name);

    template<int N>
    void InitData(string section, string field,
                  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
                  Date date, Data<T, N>& CurrentData,
                  bool interpolated = true);
    template<int N>
    void InitData(string section, string field,
                  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
                  Date date, int index, Data<T, N>& CurrentData,
                  bool interpolated = true);
    template<int N>
    void InitData(string section, string field,
                  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
                  Date date, int first_index, int second_index,
                  Data<T, N>& CurrentData);
    template<int N>
    void InitData(string input_file, Date date_min_file, T Delta_t_file,
                  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
                  Date date, Data<T, N>& CurrentData,
                  bool interpolated = true);
    template<int N>
    void InitData(string input_file, Date date_min_file, T Delta_t_file,
                  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
                  Date date, int index, Data<T, N>& CurrentData,
                  bool interpolated = true);
    template<int N>
    void InitData(string input_file, Date date_min_file, T Delta_t_file,
		  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
		  Date date, int first_index, int second_index,
		  Data<T, N>& CurrentData);
    template<int N>
    void InitData(string input_file, Date date_min_file, T Delta_t_file,
		  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
		  Date date, int index, Data<T, N>& CurrentData,
		  int Nc);
    template<int N>
    void InitData(string input_file, Date date_min_file, T Delta_t_file,
		  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
		  Date date, int first_index, int second_index, int third_index,
		  Data<T, N>& CurrentData);
    template<int N>
    void InitData(string input_file, Date date_min_file, T Delta_t_file,
		  Data<T, N>& FileData_i, Data<T, N>& FileData_f,
		  Date date, int first_index, int second_index,
		  Data<T, N>& CurrentData, int Nc);

    //LL------------------------------------------------------------
    template<int N>
    void InitData(string section, string field,
		  Data<T, N>& FileData);
    template<int N>
    void InitData(string input_file, Array<T, N>& input_data);
    //--------------------------------------------------------------
		  
    template<int N>
    void UpdateData(string section, string field, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, Data<T, N>& CurrentData,
                    bool interpolated = true);
    template<int N>
    void UpdateData(string section, string field, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, int index,
                    Data<T, N>& CurrentData, bool interpolated = true);
    template<int N>
    void UpdateData(string section, string field, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, int first_index, int second_index,
                    Data<T, N>& CurrentData);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
                    T Delta_t_file, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, Data<T, N>& CurrentData,
                    bool interpolated = true);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
                    T Delta_t_file, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, int index,
                    Data<T, N>& CurrentData, bool interpolated = true);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
		    T Delta_t_file, Data<T, N>& FileData_i,
		    Data<T, N>& FileData_f, int first_index, int second_index,
		    Data<T, N>& CurrentData);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
		    T Delta_t_file, Data<T, N>& FileData_i,
		    Data<T, N>& FileData_f, int index,
		    Data<T, N>& CurrentData, int Nc);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
		    T Delta_t_file, Data<T, N>& FileData_i,
		    Data<T, N>& FileData_f, int index,
		    Data<T, N>& CurrentData_i, Data<T, N>& CurrentData_f, int Nc);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
		    T Delta_t_file, Data<T, N>& FileData_i,
		    Data<T, N>& FileData_f, int first_index, int second_index,
		    int third_index, Data<T, N>& CurrentData_i,
		    Data<T, N>& CurrentData_f);

    template<int N>
    void UpdateData(string input_file, Date date_min_file,
		    T Delta_t_file, Data<T, N>& FileData_i,
		    Data<T, N>& FileData_f, int first_index, int second_index,
		    Data<T, N>& CurrentData, int Nc);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
		    T Delta_t_file, Data<T, N>& FileData_i,
		    Data<T, N>& FileData_f, int first_index, int second_index,
		    Data<T, N>& CurrentData_i, Data<T, N>& CurrentData_f, int Nc);
		    
    template<int N>
    void UpdateData(string section, string field, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, Data<T, N>& CurrentData_i,
                    Data<T, N>& CurrentData_f, bool interpolated = true);
    template<int N>
    void UpdateData(string section, string field, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, int index,
                    Data<T, N>& CurrentData_i, Data<T, N>& CurrentData_f,
                    bool interpolated = true);
    template<int N>
    void UpdateData(string section, string field, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, int first_index, int second_index,
                    Data<T, N>& CurrentData_i, Data<T, N>& CurrentData_f);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
                    T Delta_t_file, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, Data<T, N>& CurrentData_i,
                    Data<T, N>& CurrentData_f, bool interpolated = true);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
                    T Delta_t_file, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, int index,
                    Data<T, N>& CurrentData_i, Data<T, N>& CurrentData_f,
                    bool interpolated = true);
    template<int N>
    void UpdateData(string input_file, Date date_min_file,
                    T Delta_t_file, Data<T, N>& FileData_i,
                    Data<T, N>& FileData_f, int first_index, int second_index,
                    Data<T, N>& CurrentData_i, Data<T, N>& CurrentData_f);

    void Interpolate(Date date_min_file, T Delta_t_file, int record,
                     const Data<T, 3>& FileData_i,
                     const Data<T, 3>& FileData_f,
                     Date date, Data<T, 3>& InterpolatedData);
    void Interpolate(Date date_min_file, T Delta_t_file, int record,
                     const Data<T, 4>& FileData_i,
                     const Data<T, 4>& FileData_f,
                     Date date, Data<T, 4>& InterpolatedData);
    void Interpolate(Date date_min_file, T Delta_t_file, int record,
                     const Data<T, 2>& FileData_i,
                     const Data<T, 2>& FileData_f,
                     Date date, Data<T, 2>& InterpolatedData);
    //LL--------------------------------------------------------
    void Interpolate(Date date_min_file, T Delta_t_file, int record,
		     const Data<T, 1>& FileData_i,
		     const Data<T, 1>& FileData_f,
		     Date date, Data<T, 1>& InterpolatedData);  
    //----------------------------------------------------------
    
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_BASEMODEL_HXX
#endif
