include ../../shared.mk

all: qasm_dmsim_nvgpu 

qasm_dmsim_nvgpu: qasm_dmsim_nvgpu.cu qasm_parser.hpp ../src/dmsim_nvgpu_sin.cuh ../src/config.h 
	$(NVCC) $(NVCC_FLAGS) $(NVCC_LIBS) -D USE_NVGPU qasm_dmsim_nvgpu.cu -o $@

clean:
	rm qasm_dmsim_nvgpu

