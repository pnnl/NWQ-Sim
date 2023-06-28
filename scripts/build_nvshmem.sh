#!/bin/bash

set -e

module unload PrgEnv-nvidia
module load PrgEnv-gnu
module unload cuda
module load cudatoolkit
module load craype-accel-nvidia80

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
export FI_LOG_LEVEL=Warn

export NVSHMEM_HOME=${HOME}/nvshmem
export NVSHMEM_PREFIX=${NVSHMEM_HOME}
export NVSHMEM_USE_GDRCOPY=1
export GDRCOPY_HOME=/usr

export NVSHMEM_MPI_SUPPORT=1
export MPI_HOME=${MPICH_DIR}
export NVCC_GENCODE="-gencode=arch=compute_80,code=sm_80"

export NVSHMEM_DEFAULT_PMI2=1
export NVCUFLAGS=--allow-unsupported-compiler

export MPICC=mpicc
export CC=cc
export CXX=CC

export NVSHMEM_LIBFABRIC_SUPPORT=1
export LIBFABRIC_HOME=/opt/cray/libfabric/1.15.2.0
export LD_LIBRARY_PATH=$NVSHMEM_HOME/lib:$LD_LIBRARY_PATH

alias ls="ls --color"

export NVSHMEM_DISABLE_CUDA_VMM=1
export FI_CXI_OPTIMIZED_MRS=false
export NVSHMEM_REMOTE_TRANSPORT=libfabric

mkdir build
cd build
cmake .. -DNVSHMEM_BUILD_EXAMPLES=false -DNVSHMEM_BUILD_TESTS=false
make -j 10 install