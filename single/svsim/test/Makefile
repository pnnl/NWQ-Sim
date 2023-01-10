include ../../shared.mk

#all: qasm_svsim_benchmark_nvgpu
all: qasm_svsim_benchmark_cpu qasm_svsim_benchmark_cpu_omp qasm_svsim_benchmark_nvgpu

qasm_svsim_benchmark_nvgpu: qasm_svsim_benchmark_nvgpu.cu ../qasm/qasm_parser.hpp ../src/svsim_nvgpu_sin.cuh ../src/config.h 
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) -D USE_NVGPU qasm_svsim_benchmark_nvgpu.cu -o $@

qasm_svsim_benchmark_cpu: qasm_svsim_benchmark_cpu.cpp ../qasm/qasm_parser.hpp ../src/svsim_cpu_sin.hpp ../src/config.h 
	$(CC) $(CC_FLAGS) $(CC_LIBS) -D USE_CPU qasm_svsim_benchmark_cpu.cpp -o $@

qasm_svsim_benchmark_cpu_omp: qasm_svsim_benchmark_cpu.cpp ../qasm/qasm_parser.hpp ../src/svsim_cpu_omp.hpp ../src/config.h 
	$(CC) $(CC_FLAGS) $(CC_LIBS) -fopenmp -D USE_CPU -D USE_OMP qasm_svsim_benchmark_cpu.cpp -o $@

clean:
	rm qasm_svsim_benchmark_nvgpu qasm_svsim_benchmark_cpu qasm_svsim_benchmark_cpu_omp

