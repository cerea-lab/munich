// Copyright (C) 2009, ENPC - INRIA - EDF R&D
// Author(s): Pierre Tran
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


#ifndef POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULEPARALLEL_CXX


#include "BaseModuleParallel.hxx"

#include "BaseModule.cxx"


namespace Polyphemus
{


  //! Initializes some attributes related to the parallelization tools.
  /*! It might be always be called even when there is no parallelization
    involved. When POLYPHEMUS_PARALLEL_WITH_MPI is defined, it should be
    called after MPI_Init but before any call to other methods
    of this class. When it is not, a call to GetRank() that will return zero
    in all cases is allowed. When POLYPHEMUS_PARALLEL_WITH_OPENMP is defined,
    it should be called before any directives #pragma omp.
  */
  template<class T>
  void BaseModuleParallel::Init(BaseModel<T>& Model)
  {
#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    // Initializes MPI if not done yet.
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
      MPI_Init(NULL, NULL);

    // Initializes some attributes related to the model and the MPI universe.
    Ns_ = Model.GetNs();
    Nx_ = Model.GetNx();
    Ny_ = Model.GetNy();
    Ns_aer_ = Model.GetNs_aer();

    // Rank of the process considered in the "LAM universe".
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // Total size of the "LAM universe".
    MPI_Comm_size(MPI_COMM_WORLD, &number_slices_);
    parallelized_on_species_ = (Ns_ > 1 && Ns_ >= number_slices_);
    parallelized_on_aer_species_ = (Ns_aer_ > 1 && Ns_aer_ >= number_slices_);
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    ConfigStream config(Model.GetConfigurationFile());
    bool has_openmp_config = false;
    try
      {
        config.SetSection("[computing]");
        has_openmp_config = true;
      }
    catch (...)
      {
      }

    Nthreads_openmp_ = omp_get_max_threads();
    if (has_openmp_config)
      config.PeekValue("Number_of_threads_openmp", "integer | > 0",
                       Nthreads_openmp_);
    Nthreads_openmp_ = (Nthreads_openmp_ > omp_get_max_threads()) ?
      omp_get_max_threads() : Nthreads_openmp_;
    // The number of threads involved will remain constant.
    omp_set_dynamic(0);
#endif
  }

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
  //! Divide a 1D-interval by the number of openMP-threads.
  /*! The resulting indices comply with the FORTRAN conventions.
    \param first_index First index (under the C convention) of the interval to
    divide.
    \param last_index Last index (under the C convention) of the interval to
    divide.
    \param first_index_subseg Array of first indices (under the FORTRAN
    convention) of the resulting subsegments.
    \param last_index_subseg Array of last indices (under the FORTRAN
    convention) of the resulting subsegments.
  */
  void BaseModuleParallel::BuildSubSegment(int first_index, int last_index,
                                           Array<int, 1>& first_index_subseg,
                                           Array<int, 1>& last_index_subseg)
  {
    first_index_subseg.resize(Nthreads_openmp_);
    last_index_subseg.resize(Nthreads_openmp_);
    int dim = last_index - first_index;
    int dim_sliced = int(floor(dim / Nthreads_openmp_));
    int dim_sliced_tmp;

    for (int l = 0; l < Nthreads_openmp_; l++)
      {
        first_index_subseg(l) = (l != 0) ?
          last_index_subseg(l - 1) + 1 : first_index + 1;
        dim_sliced_tmp = (l < (dim % Nthreads_openmp_)) ?
          dim_sliced + 1 : dim_sliced;
        last_index_subseg(l) = first_index_subseg(l) - 1 + dim_sliced_tmp;
      }
  }
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
  //! Partition the domain into slices along a given direction.
  /*! The domain is sliced along the direction of dim.
    \param dim Dimension in the direction along which the partition is
    built.
    \param dim_sliced Dimension of a slice in the dim direction.
    \param first_slice_index Table of the data object starting indices of the
    slices.
    \param last_slice_index Table of the data object ending indices of the
    slices.
    \param count_slice Table of the slices extent along the sliced dimension.
    \param offset_slice Table of the buffer starting indices of the slices.
  */
  void BaseModuleParallel::BuildPartition(int dim, int& dim_sliced,
                                          Array<int, 1>& first_slice_index,
                                          Array<int, 1>& last_slice_index,
                                          Array<int, 1>& count_slice,
                                          Array<int, 1>& offset_slice)
  {
    first_slice_index.resize(number_slices_);
    last_slice_index.resize(number_slices_);
    count_slice.resize(number_slices_);
    offset_slice.resize(number_slices_);

    int l, dim_sliced_tmp;

    // Number of cells along the sliced direction.
    dim_sliced = int(floor(dim / number_slices_));
    for (l = 0; l < number_slices_; l++)
      {
        first_slice_index(l) = (l != 0) ? last_slice_index(l - 1) : 0;

        dim_sliced_tmp = (l < (dim % number_slices_)) ?
          dim_sliced + 1 : dim_sliced;

        last_slice_index(l) = first_slice_index(l) + dim_sliced_tmp;

        count_slice(l) = dim_sliced_tmp;

        offset_slice(l) = (l != 0) ?
          offset_slice(l - 1) + count_slice(l - 1) : 0;
      }

    if (rank_ < (dim % number_slices_))
      dim_sliced++;
  }


  //! Gives the edges of the current slice in a species partition.
  void BaseModuleParallel::GetEdgePartition_s(int& first_index,
                                              int& last_index)
  {
    first_index = first_slice_index_s_(rank_);
    last_index = last_slice_index_s_(rank_);
  }


  //! Gives the edges of the current slice in a aerosol species partition.
  void BaseModuleParallel::GetEdgePartition_s_aer(int& first_index,
                                                  int& last_index)
  {
    first_index = first_slice_index_s_aer_(rank_);
    last_index = last_slice_index_s_aer_(rank_);
  }


  //! Gives the edges of the current slice in a partition along x.
  void BaseModuleParallel::GetEdgePartition_x(int& first_index,
                                              int& last_index)
  {
    first_index = first_slice_index_x_(rank_);
    last_index = last_slice_index_x_(rank_);
  }


  //! Gives the edges of the current slice in a partition along y.
  void BaseModuleParallel::GetEdgePartition_y(int& first_index,
                                              int& last_index)
  {
    first_index = first_slice_index_y_(rank_);
    last_index = last_slice_index_y_(rank_);
  }


  //! Copies Array<T,3> object in an Array<T,3> permuting the axis order.
  /*! 321 gives the matching between dimensions of both objects:
    A3_out -> A3_in, 1st -> 3rd, 2nd -> 2nd, 3rd -> 1st.
  */
  template<class T>
  void BaseModuleParallel::CopyAndPermute_321(Array<T, 3>& A3_in,
                                              Array<T, 3>& A3_out)
  {
    // Loops are suited to the A3_out storage order.
    for (int i = 0; i < A3_out.extent(0); i++)
      for (int k = 0; k < A3_out.extent(1); k++)
        for (int j = 0; j < A3_out.extent(2); j++)
          A3_out(i, k, j) = A3_in(j, k, i);
  }


  //! Copies Array<T,3> object in an Array<T,3> permuting the axis order.
  /*! 213 gives the matching between dimensions of both objects:
    A3_out -> A3_in, 1st -> 2nd, 2nd -> 1st, 3rd -> 3rd.
  */
  template<class T>
  void BaseModuleParallel::CopyAndPermute_213(Array<T, 3>& A3_in,
                                              Array<T, 3>& A3_out)
  {
    // Loops are suited to the A3_out storage order.
    for (int j = 0; j < A3_out.extent(0); j++)
      for (int k = 0; k < A3_out.extent(1); k++)
        for (int i = 0; i < A3_out.extent(2); i++)
          A3_out(j, k, i) = A3_in(k, j, i);
  }


  //! Copies Array<T, 4> object in an Array<T,4> permuting the axis order.
  /*! 4231 gives the matching between dimensions of both objects:
    A4_out -> A4_in, 1st -> 4th, 2nd -> 2nd, 3rd -> 3rd, 4th -> 1st.
  */
  template<class T>
  void BaseModuleParallel::CopyAndPermute_4231(Array<T, 4>& A4_in,
                                               Array<T, 4>& A4_out)
  {
    // Loops are suited to the A4_out storage order.
    for (int i = 0; i < A4_out.extent(0); i++)
      for (int k = 0; k < A4_out.extent(1); k++)
        for (int j = 0; j < A4_out.extent(2); j++)
          for (int s = 0; s < A4_out.extent(3); s++)
            A4_out(i, k, j, s) = A4_in(s, k, j, i);
  }


  //! Copies Array<T, 5> object in an Array<T, 5> permuting the axis order.
  /*! 52341 gives the matching between dimensions of both objects:
    A5_out -> A5_in, 1st -> 5th, 2nd -> 2nd, 3rd -> 3rd, 4th -> 4th,
    5th -> 1st.
  */
  template<class T>
  void BaseModuleParallel::CopyAndPermute_52341(Array<T, 5>& A5_in,
                                                Array<T, 5>& A5_out)
  {
    // Loops are suited to the A5_out storage order.
    for (int i = 0; i < A5_out.extent(0); i++)
      for (int b = 0; b < A5_out.extent(1); b++)
        for (int k = 0; k < A5_out.extent(2); k++)
          for (int j = 0; j < A5_out.extent(3); j++)
            for (int s = 0; s < A5_out.extent(4); s++)
              A5_out(i, b, k, j, s) = A5_in(s, b, k, j, i);
  }


  //! Copies an Array<T,3> in a Array<T,3> object permuting the axis order.
  /*! 312 gives the matching between dimensions of both objects:
    A3_out -> A3_in, 1st -> 2nd, 2nd -> 3rd, 3rd -> 1st.
  */
  template<class T>
  void BaseModuleParallel::CopyAndPermute_231(Array<T, 3>& A3_in,
                                              Array<T, 3>& A3_out)
  {
    // Loops are suited to the A3 storage order.
    for (int i = 0; i < A3_in.extent(0); i++)
      for (int k = 0; k < A3_in.extent(1); k++)
        for (int j = 0; j < A3_in.extent(2); j++)
          A3_out(k, j, i) = A3_in(i, k, j);
  }


  //! Copies the array part of the current slice in an Array<T, 3> object.
  /*! Only for slave processes.  As order matters, this method should be
    used in the case of a partition of Array along its third dimension that
    is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::CopyFromSlice_x(Array<T, 3>& A3_in,
                                           Array<T, 3>& A3_out)
  {
    if (rank_ != 0)
      // Loops are suited to the A3_out storage order.
      for (int k = 0; k < A3_in.extent(1); k++)
        for (int j = 0; j < A3_in.extent(2); j++)
          for (int i = first_slice_index_x_(rank_);
               i < last_slice_index_x_(rank_); i++)
            A3_out(k, j, i) = A3_in(i, k, j);
  }


  //! Copies the array part of the current slice in an Array<T, 4> object.
  /*! Only for slave processes.  As order matters, this method should be
    used in the case of a partition of Array along its fourth dimension that
    is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::CopyFromSlice_x(Array<T, 4>& A4_in,
                                           Array<T, 4>& A4_out)
  {
    if (rank_ != 0)
      // Loops are suited to the A4_out storage order.
      for (int s = 0; s < A4_in.extent(3); s++)
        for (int k = 0; k < A4_in.extent(1); k++)
          for (int j = 0; j < A4_in.extent(2); j++)
            for (int i = first_slice_index_x_(rank_);
                 i < last_slice_index_x_(rank_); i++)
              A4_out(s, k, j, i) = A4_in(i, k, j, s);
  }


  //! Copies the array part of the current slice in an Array<T, 5> object.
  /*! Only for slave processes.  As order matters, this method should be
    used in the case of a partition of Array along its fifth dimension that
    is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::CopyFromSlice_x(Array<T, 5>& A5_in,
                                           Array<T, 5>& A5_out)
  {
    if (rank_ != 0)
      // Loops are suited to the A5_out storage order.
      for (int s = 0; s < A5_in.extent(4); s++)
        for (int b = 0; b < A5_in.extent(1); b++)
          for (int k = 0; k < A5_in.extent(2); k++)
            for (int j = 0; j < A5_in.extent(3); j++)
              for (int i = first_slice_index_x_(rank_);
                   i < last_slice_index_x_(rank_); i++)
                A5_out(s, b, k, j, i) = A5_in(i, b, k, j, s);
  }


  //! Copies the Array part of the current slice in an Array<T, 3> object.
  /*! Only for slave processes.  As order matters, this method should be
    used in the case of a partition of Array along its third dimension that
    is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::CopyToSlice_x(Array<T, 3>& A3_in,
                                         Array<T, 3>& A3_out)
  {
    // Loops are suited to the A3_in storage order.
    for (int k = 0; k < A3_out.extent(1); k++)
      for (int j = 0; j < A3_out.extent(2); j++)
        for (int i = first_slice_index_x_(rank_);
             i < last_slice_index_x_(rank_); i++)
          A3_out(i, k, j) = A3_in(k, j, i);
  }


  //! Copies the Array part of the current slice in an Array<T, 4> object.
  /*! Only for slave processes.  As order matters, this method should be
    used in the case of a partition of Array along its fourth dimension that
    is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::CopyToSlice_x(Array<T, 4>& A4_in,
                                         Array<T, 4>& A4_out)
  {
    // Loops are suited to the A4_in storage order.
    for (int s = 0; s < A4_out.extent(3); s++)
      for (int k = 0; k < A4_out.extent(1); k++)
        for (int j = 0; j < A4_out.extent(2); j++)
          for (int i = first_slice_index_x_(rank_);
               i < last_slice_index_x_(rank_); i++)
            A4_out(i, k, j, s) = A4_in(s, k, j, i);
  }


  //! Copies the Array part of the current slice in an Array<T, 5> object.
  /*! Only for slave processes.  As order matters, this method should be
    used in the case of a partition of Array along its fifth dimension that
    is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::CopyToSlice_x(Array<T, 5>& A5_in,
                                         Array<T, 5>& A5_out)
  {
    // Loops are suited to the A5_out storage order.
    for (int i = first_slice_index_x_(rank_); i < last_slice_index_x_(rank_);
         i++)
      for (int b = 0; b < A5_out.extent(1); b++)
        for (int k = 0; k < A5_out.extent(2); k++)
          for (int j = 0; j < A5_out.extent(3); j++)
            for (int s = 0; s < A5_out.extent(4); s++)
              A5_out(i, b, k, j, s) = A5_in(s, b, k, j, i);
  }


  //! Scatters an Array<T, 3> object using an Array<T, 3> object.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its third dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_x_MPI(Array<T, 3>& A3_in,
                                              Array<T, 3>& A3_out)
  {
    // Dimensions are reordered so that we can manipulate slices that
    // fit in contiguous storage locations.
    A3_out.resize(A3_in.extent(2), A3_in.extent(1), A3_in.extent(0));

    if (rank_ == 0)
      CopyAndPermute_321(A3_in, A3_out);

    MPI_Datatype TypeSlice_x;
    MPI_Type_contiguous(A3_out.extent(1) * A3_out.extent(2), MPI_DOUBLE,
                        &TypeSlice_x);
    MPI_Type_commit(&TypeSlice_x);
    MPI_Scatterv(&A3_out(0, 0, 0),
                 count_slice_x_.data(), offset_slice_x_.data(), TypeSlice_x,
                 &A3_out(first_slice_index_x_(rank_), 0, 0),
                 dim_slice_x_, TypeSlice_x, 0, MPI_COMM_WORLD);
  }


  //! Scatters an Array<T, 4> object using an Array<T, 4>.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fourth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_x_MPI(Array<T, 4>& A4_in,
                                              Array<T, 4>& A4_out)
  {
    // Dimensions are reordered so that we can manipulate slices that
    // fit in contiguous storage locations.
    A4_out.resize(A4_in.extent(3), A4_in.extent(1), A4_in.extent(2),
                  A4_in.extent(0));
    if (rank_ == 0)
      CopyAndPermute_4231(A4_in, A4_out);

    MPI_Datatype TypeSlice_x;
    MPI_Type_contiguous(A4_out.extent(1) * A4_out.extent(2)
                        * A4_out.extent(3),
                        MPI_DOUBLE, &TypeSlice_x);
    MPI_Type_commit(&TypeSlice_x);
    MPI_Scatterv(&A4_out(0, 0, 0, 0),
                 count_slice_x_.data(), offset_slice_x_.data(), TypeSlice_x,
                 &A4_out(first_slice_index_x_(rank_), 0, 0, 0),
                 dim_slice_x_, TypeSlice_x, 0, MPI_COMM_WORLD);
  }


  //! Scatters an Array<T, 5> object using an Array<T, 5>.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fifth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_x_MPI(Array<T, 5>& A5_in,
                                              Array<T, 5>& A5_out)
  {
    // Dimensions are reordered so that we can manipulate slices that
    // fit in contiguous storage locations.
    A5_out.resize(A5_in.extent(4), A5_in.extent(1), A5_in.extent(2),
                  A5_in.extent(3), A5_in.extent(0));
    if (rank_ == 0)
      CopyAndPermute_52341(A5_in, A5_out);

    MPI_Datatype TypeSlice_x;
    MPI_Type_contiguous(A5_out.extent(1) * A5_out.extent(2) * A5_out.extent(3)
                        * A5_out.extent(4),
                        MPI_DOUBLE, &TypeSlice_x);
    MPI_Type_commit(&TypeSlice_x);
    MPI_Scatterv(&A5_out(0, 0, 0, 0, 0),
                 count_slice_x_.data(), offset_slice_x_.data(), TypeSlice_x,
                 &A5_out(first_slice_index_x_(rank_), 0, 0, 0, 0),
                 dim_slice_x_, TypeSlice_x, 0, MPI_COMM_WORLD);
  }


  //! Scatters an Array<T, 3> object from master to slaves processes.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its third dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_x_MPI(Array<T, 3>& A3)
  {
    Array<T, 3> A3_tmp;
    ScatterSlice_x_MPI(A3, A3_tmp);
    CopyFromSlice_x(A3_tmp, A3);
  }


  //! Scatters an Array<T, 4> object from master to slaves processes.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fourth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_x_MPI(Array<T, 4>& A4)
  {
    Array<T, 4> A4_tmp;
    ScatterSlice_x_MPI(A4, A4_tmp);
    CopyFromSlice_x(A4_tmp, A4);
  }


  //! Scatters an Array<T, 5> object from master to slaves processes.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fifth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_x_MPI(Array<T, 5>& A5)
  {
    Array<T, 5> A5_tmp;
    ScatterSlice_x_MPI(A5, A5_tmp);
    CopyFromSlice_x(A5_tmp, A5);
  }


  //! Scatters an Array<T, 4> object from master to slaves processes.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its first dimension that is supposed to be the
    species list.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_s_MPI(Array<T, 4>& A4)
  {
    MPI_Datatype TypeSlice_s;
    MPI_Type_contiguous(A4.extent(1) * A4.extent(2) * A4.extent(3),
                        MPI_DOUBLE, &TypeSlice_s);
    MPI_Type_commit(&TypeSlice_s);
    MPI_Scatterv(&A4(0, 0, 0, 0),
                 count_slice_s_.data(), offset_slice_s_.data(), TypeSlice_s,
                 &A4(first_slice_index_s_(rank_), 0, 0, 0),
                 dim_slice_s_, TypeSlice_s, 0, MPI_COMM_WORLD);
  }


  //! Scatters an Array<T, 5> object from master to slaves processes.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its first dimension that is supposed to be the
    aerosol species list.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_s_aer_MPI(Array<T, 5>& A5)
  {
    MPI_Datatype TypeSlice_s_aer;
    MPI_Type_contiguous(A5.extent(1) * A5.extent(2) * A5.extent(3)
                        * A5.extent(4),
                        MPI_DOUBLE, &TypeSlice_s_aer);
    MPI_Type_commit(&TypeSlice_s_aer);
    MPI_Scatterv(&A5(0, 0, 0, 0, 0),
                 count_slice_s_aer_.data(), offset_slice_s_aer_.data(),
                 TypeSlice_s_aer,
                 &A5(first_slice_index_s_aer_(rank_), 0, 0, 0, 0),
                 dim_slice_s_aer_, TypeSlice_s_aer, 0, MPI_COMM_WORLD);
  }


  //! Scatters an Array<T, 3> object using an Array<T, 3>.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its second dimension that is supposed to be Y.
  */
  template<class T>
  void BaseModuleParallel::ScatterSlice_y_MPI(Array<T, 3>& A3_in,
                                              Array<T, 3>& A3_out)
  {
    // Dimensions are reordered so that we can manipulate slices that
    // fit in contiguous storage locations.
    A3_out.resize(A3_in.extent(1), A3_in.extent(0), A3_in.extent(2));

    if (rank_ == 0)
      CopyAndPermute_213(A3_in, A3_out);

    MPI_Datatype TypeSlice_y;
    MPI_Type_contiguous(A3_out.extent(1) * A3_out.extent(2),
                        MPI_DOUBLE, &TypeSlice_y);
    MPI_Type_commit(&TypeSlice_y);
    MPI_Scatterv(&A3_out(0, 0, 0),
                 count_slice_y_.data(), offset_slice_y_.data(), TypeSlice_y,
                 &A3_out(first_slice_index_y_(rank_), 0, 0),
                 dim_slice_y_, TypeSlice_y, 0, MPI_COMM_WORLD);
  }


  //! Gathers a Array<T, 3> object using an Array<T, 3>.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its third dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_x_MPI(Array<T, 3>& A3_in,
                                             Array<T, 3>& A3_out)
  {
    MPI_Datatype TypeSlice_x;
    MPI_Type_contiguous(A3_in.extent(1) * A3_in.extent(2),
                        MPI_DOUBLE, &TypeSlice_x);
    MPI_Type_commit(&TypeSlice_x);
    MPI_Gatherv(&A3_in(first_slice_index_x_(rank_), 0, 0),
                dim_slice_x_, TypeSlice_x,
                &A3_in(0, 0, 0),
                count_slice_x_.data(), offset_slice_x_.data(), TypeSlice_x,
                0, MPI_COMM_WORLD);

    if (rank_ == 0)
      CopyAndPermute_321(A3_in, A3_out);
  }


  //! Gathers an Array<T, 4> object using an Array<T, 4>.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fourth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_x_MPI(Array<T, 4>& A4_in,
                                             Array<T, 4>& A4_out)
  {
    MPI_Datatype TypeSlice_x;
    MPI_Type_contiguous(A4_in.extent(1) * A4_in.extent(2) * A4_in.extent(3),
                        MPI_DOUBLE, &TypeSlice_x);
    MPI_Type_commit(&TypeSlice_x);
    MPI_Gatherv(&A4_in(first_slice_index_x_(rank_), 0, 0, 0),
                dim_slice_x_, TypeSlice_x,
                &A4_in(0, 0, 0, 0),
                count_slice_x_.data(), offset_slice_x_.data(), TypeSlice_x,
                0, MPI_COMM_WORLD);

    if (rank_ == 0)
      CopyAndPermute_4231(A4_in, A4_out);
  }


  //! Gathers an Array<T, 5> object using an Array<T, 5>.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fifth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_x_MPI(Array<T, 5>& A5_in,
                                             Array<T, 5>& A5_out)
  {
    MPI_Datatype TypeSlice_x;
    MPI_Type_contiguous(A5_in.extent(1) * A5_in.extent(2) * A5_in.extent(3)
                        * A5_in.extent(4),
                        MPI_DOUBLE, &TypeSlice_x);
    MPI_Type_commit(&TypeSlice_x);
    MPI_Gatherv(&A5_in(first_slice_index_x_(rank_), 0, 0, 0, 0),
                dim_slice_x_, TypeSlice_x,
                &A5_in(0, 0, 0, 0, 0),
                count_slice_x_.data(), offset_slice_x_.data(), TypeSlice_x,
                0, MPI_COMM_WORLD);

    if (rank_ == 0)
      CopyAndPermute_52341(A5_in, A5_out);
  }


  //! Gathers an Array<T, 3> object from slaves to master process.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its third dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_x_MPI(Array<T, 3>& A3)
  {
    Array<T, 3> A3_tmp(A3.extent(2), A3.extent(0), A3.extent(1));
    CopyToSlice_x(A3, A3_tmp);
    GatherSlice_x_MPI(A3_tmp, A3);
  }


  //! Gathers an Array<T, 4> object from slaves to master process.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fourth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_x_MPI(Array<T, 4>& A4)
  {
    Array<T, 4> A4_tmp(A4.extent(3), A4.extent(1), A4.extent(2),
                       A4.extent(0));
    CopyToSlice_x(A4, A4_tmp);
    GatherSlice_x_MPI(A4_tmp, A4);
  }


  //! Gathers an Array<T, 5> object from slaves to master process.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its fifth dimension that is supposed to be X.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_x_MPI(Array<T, 5>& A5)
  {
    Array<T, 5> A5_tmp(A5.extent(4), A5.extent(1), A5.extent(2),
                       A5.extent(3), A5.extent(0));
    CopyToSlice_x(A5, A5_tmp);
    GatherSlice_x_MPI(A5_tmp, A5);
  }


  //! Gathers an Array<T, 4> object from slaves to master process.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its first dimension that is supposed to be the
    species list.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_s_MPI(Array<T, 4>& A4)
  {
    MPI_Datatype TypeSlice_s;
    MPI_Type_contiguous(A4.extent(1) * A4.extent(2) * A4.extent(3),
                        MPI_DOUBLE, &TypeSlice_s);
    MPI_Type_commit(&TypeSlice_s);
    MPI_Gatherv(&A4(first_slice_index_s_(rank_), 0, 0, 0),
                dim_slice_s_, TypeSlice_s,
                &A4(0, 0, 0, 0),
                count_slice_s_.data(), offset_slice_s_.data(), TypeSlice_s,
                0, MPI_COMM_WORLD);
  }


  //! Gathers an Array<T, 5> object from slaves to master process.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its first dimension that is supposed to be the
    aerosol species list.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_s_aer_MPI(Array<T, 5>& A5)
  {
    MPI_Datatype TypeSlice_s_aer;
    MPI_Type_contiguous(A5.extent(1) * A5.extent(2) * A5.extent(3)
                        * A5.extent(4),
                        MPI_DOUBLE, &TypeSlice_s_aer);
    MPI_Type_commit(&TypeSlice_s_aer);
    MPI_Gatherv(&A5(first_slice_index_s_aer_(rank_), 0, 0, 0, 0),
                dim_slice_s_aer_, TypeSlice_s_aer,
                &A5(0, 0, 0, 0, 0),
                count_slice_s_aer_.data(), offset_slice_s_aer_.data(),
                TypeSlice_s_aer,
                0, MPI_COMM_WORLD);
  }


  //! Gathers an Array<T, 3> object using an Array<T, 3>.
  /*! As order matters, this method should be used in the case of a
    partition of Array along its second dimension that is supposed to be Y.
  */
  template<class T>
  void BaseModuleParallel::GatherSlice_y_MPI(Array<T, 3>& A3_in,
                                             Array<T, 3>& A3_out)

  {
    MPI_Datatype TypeSlice_y;
    MPI_Type_contiguous(A3_in.extent(1) * A3_in.extent(2),
                        MPI_DOUBLE, &TypeSlice_y);
    MPI_Type_commit(&TypeSlice_y);
    MPI_Gatherv(&A3_in(first_slice_index_y_(rank_), 0, 0),
                dim_slice_y_, TypeSlice_y,
                &A3_in(0, 0, 0),
                count_slice_y_.data(), offset_slice_y_.data(), TypeSlice_y,
                0, MPI_COMM_WORLD);

    if (rank_ == 0)
      CopyAndPermute_213(A3_in, A3_out);
  }
#endif


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULEPARALLEL_CXX
#endif
