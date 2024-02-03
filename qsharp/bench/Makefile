include ./shared.mk

all: Grover

Grover: Grover_driver.cc Grover.ll ../qsharp_wrapper.cu ../dump_qasm.hpp
	$(QIRCC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) -o Grover.o -c Grover.ll
	$(CC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) $(SCALE) -D USE_NVGPU -o Grover_driver.o -c Grover_driver.cc
	$(NVCC) -ccbin $(CC) $(NVCC_FLAGS) $(QIR_BRIDGE_FLAGS) $(NVCC_LIBS) -D USE_NVGPU ../qsharp_wrapper.cu Grover.o Grover_driver.o -o $@

clean:
	rm *.o Grover

