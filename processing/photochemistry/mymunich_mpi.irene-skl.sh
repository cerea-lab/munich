#!/bin/bash

export LD_LIBRARY_PATH=''
module purge
module switch dfldatadir/shr10274
 module load flavor/buildcompiler/gcc/4.8
#module load intel/17.0.4.196
#module load flavor/hdf5/parallel
module load mpi/openmpi
#module load hdf5
#module load netcdf-fortran
#module load grib
#module load jasper
#module load blitz
module load blitz/1.0.2
module load scons

#source ${FORTRAN_INTEL_ROOT}/bin/compilervars.sh intel64
#source ${C_INTEL_ROOT}/bin/compilervars.sh intel64


#export LD_LIBRARY_PATH=/ccc/products/blitz-0.10/intel--20.0.0__openmpi--4.0.1/default/include:$LD_LIBRARY_PATH
#export CPLUS_INCLUDE_PATH=/ccc/products/blitz-0.10/intel--20.0.0__openmpi--4.0.1/default/include
export CPLUS_INCLUDE_PATH=/ccc/products/blitz-1.0.2/system/default/include:$CPLUS_INCLUDE_PATH
#export CPLUS_INCLUDE_PATH=/ccc/products/blitz-1.0.2/intel--17.0.4.196/default/include
#export LD_LIBRARY_PATH=/usr/lib/:$LD_LIBRARY_PATH
