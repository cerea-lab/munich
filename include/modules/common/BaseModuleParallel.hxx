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


#ifndef POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULEPARALLEL_HXX

#include "BaseModule.hxx"

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
#include <omp.h>
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
#include <mpi.h>
#endif


namespace Polyphemus
{


  ////////////////////////
  // BASEMODULEPARALLEL //
  ////////////////////////


  //! This class provides tools for parallelization with MPI.
  /*! Modules to parallelize should be derived from this class.
   */
  class BaseModuleParallel: public BaseModule
  {


  private:
#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int Nthreads_openmp_;
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    int Nx_;
    int Ny_;
    int Ns_;
    int Ns_aer_;

    int rank_;
    int number_slices_;
    bool parallelized_on_species_;
    bool parallelized_on_aer_species_;

    Array<int, 1> first_slice_index_x_;
    Array<int, 1> last_slice_index_x_;
    Array<int, 1> count_slice_x_;
    Array<int, 1> offset_slice_x_;
    int dim_slice_x_;

    Array<int, 1> first_slice_index_y_;
    Array<int, 1> last_slice_index_y_;
    Array<int, 1> count_slice_y_;
    Array<int, 1> offset_slice_y_;
    int dim_slice_y_;

    Array<int, 1> first_slice_index_s_;
    Array<int, 1> last_slice_index_s_;
    Array<int, 1> count_slice_s_;
    Array<int, 1> offset_slice_s_;
    int dim_slice_s_;

    Array<int, 1> first_slice_index_s_aer_;
    Array<int, 1> last_slice_index_s_aer_;
    Array<int, 1> count_slice_s_aer_;
    Array<int, 1> offset_slice_s_aer_;
    int dim_slice_s_aer_;

    void BuildPartition(int dim, int& dim_para,
                        Array<int, 1>& first, Array<int, 1>& last,
                        Array<int, 1>& count, Array<int, 1>& offset);

    template<class T>
    void CopyAndPermute_321(Array<T, 3>& A3_in, Array<T, 3>& A3_out);
    template<class T>
    void CopyAndPermute_213(Array<T, 3>& A3_in, Array<T, 3>& A3_out);
    template<class T>
    void CopyAndPermute_4231(Array<T, 4>& A4_in, Array<T, 4>& A4_out);
    template<class T>
    void CopyAndPermute_52341(Array<T, 5>& A5_in, Array<T, 5>& A5_out);
    template<class T>
    void CopyAndPermute_231(Array<T, 3>& A3_in, Array<T, 3>& A3_out);
#endif

  public:
    template<class T>
    void Init(BaseModel<T>& Model);

#ifdef POLYPHEMUS_PARALLEL_WITH_OPENMP
    int GetNthreads_openmp()
    {
      return Nthreads_openmp_;
    }
    void BuildSubSegment(int first_index, int last_index,
                         Array<int, 1>& first_index_subslice,
                         Array<int, 1>& last_index_subslice);
#endif

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    int GetRank()
    {
      return rank_;
    }
    bool ParallelizedOnSpecies()
    {
      return parallelized_on_species_;
    }
    bool ParallelizedOnAerSpecies()
    {
      return parallelized_on_aer_species_;
    }

    void BuildPartition_s()
    {
      BuildPartition(Ns_, dim_slice_s_,
                     first_slice_index_s_, last_slice_index_s_,
                     count_slice_s_, offset_slice_s_);
    }

    void BuildPartition_s_aer()
    {
      BuildPartition(Ns_aer_, dim_slice_s_aer_,
                     first_slice_index_s_aer_, last_slice_index_s_aer_,
                     count_slice_s_aer_, offset_slice_s_aer_);
    }

    void BuildPartition_x()
    {
      BuildPartition(Nx_, dim_slice_x_,
                     first_slice_index_x_, last_slice_index_x_,
                     count_slice_x_, offset_slice_x_);
    }

    void BuildPartition_y()
    {
      BuildPartition(Ny_, dim_slice_y_,
                     first_slice_index_y_, last_slice_index_y_,
                     count_slice_y_, offset_slice_y_);
    }

    void GetEdgePartition_s(int& first_index, int& last_index);
    void GetEdgePartition_s_aer(int& first_index, int& last_index);
    void GetEdgePartition_x(int& first_index, int& last_index);
    void GetEdgePartition_y(int& first_index, int& last_index);

    template<class T>
    void CopyFromSlice_x(Array<T, 3>& A3_int, Array<T, 3>& A3_out);
    template<class T>
    void CopyFromSlice_x(Array<T, 3>& A3, Data<T, 3>& D3)
    {
      CopyFromSlice_x(A3, D3.GetArray());
    }

    template<class T>
    void CopyFromSlice_x(Array<T, 4>& A4_int, Array<T, 4>& A4_out);
    template<class T>
    void CopyFromSlice_x(Array<T, 4>& A4, Data<T, 4>& D4)
    {
      CopyFromSlice_x(A4, D4.GetArray());
    }

    template<class T>
    void CopyFromSlice_x(Array<T, 5>& A5_int, Array<T, 5>& A5_out);
    template<class T>
    void CopyFromSlice_x(Array<T, 5>& A5, Data<T, 5>& D5)
    {
      CopyFromSlice_x(A5, D5.GetArray());
    }

    template<class T>
    void CopyToSlice_x(Array<T, 3>& A3_int, Array<T, 3>& A3_out);
    template<class T>
    void CopyToSlice_x(Data<T, 3>& D3, Array<T, 3>& A3)
    {
      CopyToSlice_x(D3.GetArray(), A3);
    }

    template<class T>
    void CopyToSlice_x(Array<T, 4>& A4_int, Array<T, 4>& A4_out);
    template<class T>
    void CopyToSlice_x(Data<T, 4>& D4, Array<T, 4>& A4)
    {
      CopyToSlice_x(D4.GetArray(), A4);
    }

    template<class T>
    void CopyToSlice_x(Array<T, 5>& A5_int, Array<T, 5>& A5_out);
    template<class T>
    void CopyToSlice_x(Data<T, 5>& D5, Array<T, 5>& A5)
    {
      CopyToSlice_x(D5.GetArray(), A5);
    }

    template<class T>
    void ScatterSlice_x_MPI(Array<T, 3>& A3_int, Array<T, 3>& A3_out);
    template<class T>
    void ScatterSlice_x_MPI(Data<T, 3>& D3, Array<T, 3>& A3)
    {
      ScatterSlice_x_MPI(D3.GetArray(), A3);
    }

    template<class T>
    void ScatterSlice_x_MPI(Array<T, 4>& A4_int, Array<T, 4>& A4_out);
    template<class T>
    void ScatterSlice_x_MPI(Data<T, 4>& D4, Array<T, 4>& A4)
    {
      ScatterSlice_x_MPI(D4.GetArray(), A4);
    }

    template<class T>
    void ScatterSlice_x_MPI(Array<T, 5>& A5_int, Array<T, 5>& A5_out);
    template<class T>
    void ScatterSlice_x_MPI(Data<T, 5>& D5, Array<T, 5>& A5)
    {
      ScatterSlice_x_MPI(D5.GetArray(), A5);
    }

    template<class T>
    void ScatterSlice_x_MPI(Array<T, 3>& A3);
    template<class T>
    void ScatterSlice_x_MPI(Data<T, 3>& D3)
    {
      ScatterSlice_x_MPI(D3.GetArray());
    }

    template<class T>
    void ScatterSlice_x_MPI(Array<T, 4>& A4);
    template<class T>
    void ScatterSlice_x_MPI(Data<T, 4>& D4)
    {
      ScatterSlice_x_MPI(D4.GetArray());
    }


    template<class T>
    void ScatterSlice_x_MPI(Array<T, 5>& A5);
    template<class T>
    void ScatterSlice_x_MPI(Data<T, 5>& D5)
    {
      ScatterSlice_x_MPI(D5.GetArray());
    }

    template<class T>
    void ScatterSlice_s_MPI(Array<T, 4>& A4);
    template<class T>
    void ScatterSlice_s_MPI(Data<T, 4>& D4)
    {
      ScatterSlice_s_MPI(D4.GetArray());
    }

    template<class T>
    void ScatterSlice_s_aer_MPI(Array<T, 5>& A5);
    template<class T>
    void ScatterSlice_s_aer_MPI(Data<T, 5>& D5)
    {
      ScatterSlice_s_aer_MPI(D5.GetArray());
    }

    template<class T>
    void ScatterSlice_y_MPI(Array<T, 3>& A3_int, Array<T, 3>& A3_out);
    template<class T>
    void ScatterSlice_y_MPI(Data<T, 3>& D3, Array<T, 3>& A3)
    {
      ScatterSlice_y_MPI(D3.GetArray(), A3);
    }

    template<class T>
    void GatherSlice_x_MPI(Array<T, 3>& A3_in, Array<T, 3>& A3_out);
    template<class T>
    void GatherSlice_x_MPI(Array<T, 3>& A3, Data<T, 3>& D3)
    {
      GatherSlice_x_MPI(A3, D3.GetArray());
    }

    template<class T>
    void GatherSlice_x_MPI(Array<T, 4>& A4_in, Array<T, 4>& A4_out);
    template<class T>
    void GatherSlice_x_MPI(Array<T, 4>& A4, Data<T, 4>& D4)
    {
      GatherSlice_x_MPI(A4, D4.GetArray());
    }

    template<class T>
    void GatherSlice_x_MPI(Array<T, 5>& A5_in, Array<T, 5>& A5_out);
    template<class T>
    void GatherSlice_x_MPI(Array<T, 5>& A5, Data<T, 5>& D5)
    {
      GatherSlice_x_MPI(A5, D5.GetArray());
    }

    template<class T>
    void GatherSlice_x_MPI(Array<T, 3>& A3);
    template<class T>
    void GatherSlice_x_MPI(Data<T, 3>& D3)
    {
      GatherSlice_x_MPI(D3.GetArray());
    }

    template<class T>
    void GatherSlice_x_MPI(Array<T, 4>& A4);
    template<class T>
    void GatherSlice_x_MPI(Data<T, 4>& D4)
    {
      GatherSlice_x_MPI(D4.GetArray());
    }

    template<class T>
    void GatherSlice_x_MPI(Array<T, 5>& A5);
    template<class T>
    void GatherSlice_x_MPI(Data<T, 5>& D5)
    {
      GatherSlice_x_MPI(D5.GetArray());
    }

    template<class T>
    void GatherSlice_s_MPI(Array<T, 4>& A4);
    template<class T>
    void GatherSlice_s_MPI(Data<T, 4>& D4)
    {
      GatherSlice_s_MPI(D4.GetArray());
    }

    template<class T>
    void GatherSlice_s_aer_MPI(Array<T, 5>& A5);
    template<class T>
    void GatherSlice_s_aer_MPI(Data<T, 5>& D5)
    {
      GatherSlice_s_aer_MPI(D5.GetArray());
    }

    template<class T>
    void GatherSlice_y_MPI(Array<T, 3>& A3_in, Array<T, 3>& A3_out);
    template<class T>
    void GatherSlice_y_MPI(Array<T, 3>& A3, Data<T, 3>& D3)
    {
      GatherSlice_y_MPI(A3, D3.GetArray());
    }
#else
    int GetRank()
    {
      return 0;
    }
#endif
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_COMMON_BASEMODULEPARALLEL_HXX
#endif
