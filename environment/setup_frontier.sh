#Environment settings for ORNL Frontier HPC
module load cmake
module load rocm/6.4.2
module load cray-mpich
module load craype-accel-amd-gfx90a
module load python/3.13.0

export cc=gcc
export CC=g++

export MY_HIP_ARCH=gfx90a
export NWQSIM_ENABLE_HIP=ON
export MPICH_GPU_SUPPORT_ENABLED=0
