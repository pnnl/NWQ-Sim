#Environment settings for TACC H100

module load cmake
module load cuda/12.0

export MY_CUDA_ARCH=90
export cc=gcc
export CC=g++

