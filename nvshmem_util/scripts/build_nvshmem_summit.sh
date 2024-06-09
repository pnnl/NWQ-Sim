#!/bin/bash

set -e

module load cmake
module load gcc/11.4.0
module load cuda/12.2.0
module load spectrum-mpi

export NVSHMEM_HOME=${HOME}/nvshmem
export NVSHMEM_PREFIX=${NVSHMEM_HOME}
export NVSHMEM_USE_GDRCOPY=0

export NVSHMEM_MPI_SUPPORT=1
export MPI_HOME=${MPI_ROOT}
export NVCC_GENCODE="-gencode=arch=compute_70,code=sm_70"

#export NVSHMEM_DEFAULT_PMI2=1
export NVCUFLAGS=--allow-unsupported-compiler

export MPICC=mpicc
export CC=cc
export CXX=c++

export LD_LIBRARY_PATH=$NVSHMEM_HOME/lib:$LD_LIBRARY_PATH

alias ls="ls --color"

export NVSHMEM_DISABLE_CUDA_VMM=1
export FI_CXI_OPTIMIZED_MRS=false

mkdir build
cd build
cmake .. -DNVSHMEM_BUILD_EXAMPLES=false -DNVSHMEM_BUILD_TESTS=false -DCUDA_HOME=$CUDA_DIR
make -j 10 install
