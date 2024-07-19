#Environment settings for ORNL Frontier HPC
module load cmake
module load rocm/5.4.0
module load cray-mpich
module load craype-accel-amd-gfx90a

export cc=gcc
export CC=g++

export MY_HIP_ARCH=gfx90a
export MPICH_GPU_SUPPORT_ENABLED=1
