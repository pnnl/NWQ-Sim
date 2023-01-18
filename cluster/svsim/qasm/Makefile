include ../../shared.mk

all: qasm_svsim_nvgpu
#all: qasm_svsim_cpu

qasm_svsim_nvgpu: qasm_svsim_nvgpu.cu qasm_parser.hpp ../src/svsim_nvgpu_mpi.cuh ../src/config.h 
	$(NVCC) -ccbin mpicxx -lnvshmem $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) -D USE_NVGPU qasm_svsim_nvgpu.cu -o $@

qasm_svsim_cpu: qasm_svsim_cpu.cpp qasm_parser.hpp ../src/svsim_cpu_mpi.hpp ../src/config.h
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(SCALE) -D USE_CPU qasm_svsim_cpu.cpp -o $@

clean:
	rm qasm_svsim_nvgpu qasm_svsim_cpu

