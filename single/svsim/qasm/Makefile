include ../../shared.mk

all: qasm_svsim_nvgpu 
#all: qasm_svsim_cpu_omp qasm_svsim_nvgpu qasm_svsim_cpu_avx512 qasm_svsim_cpu qasm_svsim_amdgpu

qasm_svsim_nvgpu: qasm_svsim_nvgpu.cu qasm_parser.hpp ../src/svsim_nvgpu_sin.cuh ../src/config.h 
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LIBS) -D USE_NVGPU qasm_svsim_nvgpu.cu -o $@

qasm_svsim_cpu: qasm_svsim_cpu.cpp qasm_parser.hpp ../src/svsim_cpu_sin.hpp ../src/config.h 
	$(CC) $(CC_FLAGS) $(CC_LIBS) -D USE_CPU qasm_svsim_cpu.cpp -o $@

qasm_svsim_cpu_avx512: qasm_svsim_cpu.cpp qasm_parser.hpp ../src/svsim_cpu_sin.hpp ../src/config.h  ../src/svsim_cpu_avx512.hpp
	$(CC) $(CC_FLAGS) $(CC_LIBS) -D USE_CPU -D USE_AVX512 -mavx512dq -mavx512vl qasm_svsim_cpu.cpp -o $@

qasm_svsim_cpu_omp: qasm_svsim_cpu.cpp qasm_parser.hpp ../src/svsim_cpu_omp.hpp ../src/config.h 
	$(CC) $(CC_FLAGS) $(CC_LIBS) -fopenmp -D USE_CPU -D USE_OMP qasm_svsim_cpu.cpp -o $@

qasm_svsim_amdgpu: qasm_svsim_amdgpu.cpp qasm_parser.hpp ../src/svsim_amdgpu_sin.hpp ../src/config.h 
	$(HIPCC) $(HIPCC_FLAGS) $(HIPCC_LIBS) -D USE_AMDGPU qasm_svsim_amdgpu.cpp -o $@

clean:
	rm qasm_svsim_cpu_omp qasm_svsim_nvgpu qasm_svsim_cpu_avx512 qasm_svsim_cpu qasm_svsim_amdgpu

