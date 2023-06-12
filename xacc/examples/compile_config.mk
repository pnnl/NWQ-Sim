

###################### NVIDIA GPU Configuration #######################
######################         CLUSTER          #######################
NVCC_CLUSTER = /opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/bin/nvcc
NVCC_FLAGS_CLUSTER = -O3 -arch=sm_80 -m64 -allow-unsupported-compiler -std=c++17 -rdc=true --compiler-options -fPIC -ccbin CC
NVCC_LIBS_CLUSTER = -lm -lcuda -lfabric -I$(HOME)/nvshmem/include -L$(HOME)/nvshmem/lib/ -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64 -L/opt/cray/libfabric/1.15.2.0/lib64/ -lnvidia-ml -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64/stubs/
######################         SINGLE           #######################
NVCC_SINGLE = /opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/bin/nvcc
NVCC_FLAGS_SINGLE = -O3 -arch=sm_80 -m64 -std=c++17 -rdc=true --compiler-options -fPIC
NVCC_LIBS_SINGLE = 
#######################################################################

####################### Host CPU Configuration ########################
CC = g++
CC_FLAGS = -O3 -m64 -std=c++17 -fPIC -fopenmp -w
CC_LIBS = 
#######################################################################
