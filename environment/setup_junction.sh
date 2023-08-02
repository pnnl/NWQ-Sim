module load rocm/5.3.0 
module load gcc/10.3.0

export cc=gcc
export CC=g++

export MY_HIP_ARCH=gfx908
export MY_CUDA_ARCH=70