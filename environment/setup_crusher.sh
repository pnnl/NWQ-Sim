#Environment settings for ORNL Crusher HPC
module load cmake
module load rocm/5.3.0
module load cray-mpich
module load craype-accel-amd-gfx90a

export cc=gcc
export CC=g++

export MPICH_GPU_SUPPORT_ENABLED=1
export MY_HIP_ARCH=gfx90a