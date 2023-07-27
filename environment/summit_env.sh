#Environment settings for ORNL Summit HPC

module load cmake
module load gcc/9.3.0
module load cuda/11.5.2
module load gdrcopy/2.3
module load spectrum-mpi

export SHMEM_SYMMETRIC_HEAP_SIZE=2G
export NVSHMEM_SYMMETRIC_SIZE=8G

export LD_LIBRARY_PATH=$HOME/nvshmem/lib:$LD_LIBRARY_PATH

export cc=gcc
export CC=g++

export MY_CUDA_ARCH=70