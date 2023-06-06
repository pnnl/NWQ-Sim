###################### NVIDIA GPU Configuration #######################
NVCC = /opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/bin/nvcc
NVCC_FLAGS = -O3 --generate-code=arch=compute_80,code=[compute_80,sm_80] -std=c++17 -rdc=true --compiler-options -fPIC -ccbin mpicxx 
NVCC_LIBS = -lm -lcuda -lfabric -lnvidia-ml -lnvshmem -I/$(HOME)/nvshmem/include -L/$(HOME)/nvshmem/lib/ -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64 -L/opt/cray/libfabric/1.15.2.0/lib64/ -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64/stubs/


nwq_qasm_nvgpu_mpi: nwq_qasm_nvgpu_mpi.cu
	$(NVCC)  $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) nwq_qasm_nvgpu_mpi.cu -o $@

clean:
	rm -f nwq_qasm_nvgpu_mpi