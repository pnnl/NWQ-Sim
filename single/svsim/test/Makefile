include ../../shared.mk

all: qasm_svsim_benchmark_nvgpu

qasm_svsim_benchmark_nvgpu: qasm_svsim_benchmark_nvgpu.cu ../qasm/qasm_parser.hpp ../src/svsim_nvgpu_sin.cuh ../src/config.h 
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LIBS) $(SCALE) -D USE_NVGPU qasm_svsim_benchmark_nvgpu.cu -o $@

clean:
	rm qasm_svsim_benchmark_nvgpu

