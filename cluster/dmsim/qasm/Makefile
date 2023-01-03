include ../../shared.mk

all: qasm_dmsim_nvgpu

qasm_dmsim_nvgpu: qasm_dmsim_nvgpu.cu qasm_parser.hpp ../src/dmsim_nvgpu_mpi.cuh ../src/config.h 
	$(NVCC) -ccbin mpicxx -lnvshmem $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) -D USE_NVGPU qasm_dmsim_nvgpu.cu -o $@

clean:
	rm qasm_dmsim_nvgpu

