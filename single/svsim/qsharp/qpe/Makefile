include ../../../shared.mk

all: qpe_nvgpu

qpe_nvgpu: main.cpp est-energy.o JW_QIR.o JW_json.o main.o ../../src/svsim_nvgpu_sin.cuh ../../src/config.h ../qsharp_wrapper.cu 
	$(NVCC) -ccbin $(CC) $(NVCC_FLAGS) $(QIR_BRIDGE_FLAGS) $(NVCC_LIBS) -D USE_NVGPU ../qsharp_wrapper.cu est-energy.o main.o JW_QIR.o JW_json.o -o $@

est-energy.o: est-energy.IBMQ.ll
	$(QIRCC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) -o est-energy.o -c est-energy.IBMQ.ll

main.o: main.cpp
	$(CC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) $(SCALE) -D USE_NVGPU -o main.o -c main.cpp

JW_QIR.o: JW_QIR.cpp
	$(QIRCC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) $(SCALE) -D USE_NVGPU -o JW_QIR.o -c JW_QIR.cpp

JW_json.o: JW_json.cpp
	$(CC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) $(SCALE) -D USE_NVGPU -o JW_json.o -c JW_json.cpp

clean:
	rm *.o qpe_nvgpu
