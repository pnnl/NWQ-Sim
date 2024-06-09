#Environment settings for ORNL Summit HPC

module load cmake
module load gcc/11.4.0
module load cuda/12.2.0
module load spectrum-mpi
module load openblas

export SHMEM_SYMMETRIC_HEAP_SIZE=2G
export NVSHMEM_SYMMETRIC_SIZE=8G

export LD_LIBRARY_PATH=$HOME/nvshmem/lib:$LD_LIBRARY_PATH

export cc=gcc
export CC=g++

export MY_CUDA_ARCH=70
