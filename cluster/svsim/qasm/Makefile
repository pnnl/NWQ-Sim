include ../../shared.mk

all: qasm_svsim_nvgpu

qasm_svsim_nvgpu: qasm_svsim_nvgpu.cu qasm_parser.hpp ../src/svsim_nvgpu_mpi.cuh ../src/config.h 
	$(NVCC) -ccbin mpicxx -lnvshmem $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) -D USE_NVGPU qasm_svsim_nvgpu.cu -o $@

clean:
	rm qasm_svsim_nvgpu

