include ../../shared.mk

all: qasm_svsim_nvgpu 

qasm_svsim_nvgpu: qasm_svsim_nvgpu.cu qasm_parser.hpp ../src/svsim_nvgpu_sin.cuh ../src/config.h 
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LIBS) -D USE_NVGPU qasm_svsim_nvgpu.cu -o $@

clean:
	rm qasm_svsim_nvgpu

