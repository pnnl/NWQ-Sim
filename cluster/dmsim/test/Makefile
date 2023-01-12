include ../../shared.mk

#all: qasm_dmsim_benchmark_nvgpu
all: qasm_dmsim_benchmark_cpu


qasm_dmsim_benchmark_nvgpu: qasm_dmsim_benchmark_nvgpu.cu ../qasm/qasm_parser.hpp ../src/dmsim_nvgpu_mpi.cuh ../src/config.h 
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) -D USE_NVGPU qasm_dmsim_benchmark_nvgpu.cu -o $@

qasm_dmsim_benchmark_cpu: qasm_dmsim_benchmark_cpu.cpp ../qasm/qasm_parser.hpp ../src/dmsim_cpu_mpi.hpp ../src/config.h 
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(SCALE) -D USE_CPU qasm_dmsim_benchmark_cpu.cpp -o $@

clean:
	rm qasm_dmsim_benchmark_nvgpu qasm_dmsim_benchmark_cpu

