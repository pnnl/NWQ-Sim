###################### NVIDIA GPU Configuration #######################
NVCC = /opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/bin/nvcc
NVCC_FLAGS = -O3 -arch=sm_80 -m64 -allow-unsupported-compiler -std=c++17 -rdc=true --compiler-options -fPIC -ccbin CC
NVCC_LIBS = -lm -lcuda -lfabric -I/$(HOME)/nvshmem/include -L/$(HOME)/nvshmem/lib/ -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64 -L/opt/cray/libfabric/1.15.2.0/lib64/ -lnvidia-ml -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64/stubs/
#######################################################################

####################### Host CPU Configuration ########################
CC = CC
CC_FLAGS = -O3 -m64 -std=c++14 -fPIC -fopenmp
CC_LIBS = 
#######################################################################

all: nwq_qasm_nvgpu_mpi

nwq_qasm_nvgpu_mpi: nwq_qasm_nvgpu_mpi.cu
	$(NVCC) -ccbin mpicxx -lnvshmem $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) nwq_qasm_nvgpu_mpi.cu -o $@

clean:
	rm -f nwq_qasm_nvgpu_mpi