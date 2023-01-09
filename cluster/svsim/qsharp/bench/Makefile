include ../../../shared.mk

all: TeleportationSample_nvgpu Qrng_nvgpu Grover_nvgpu

common: ../../src/svsim_nvgpu_mpi.cuh ../../src/config.h ../qsharp_wrapper.cu

Grover_nvgpu: Grover_driver.cc Grover.ll $(common)
	$(QIRCC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) -o Grover.o -c Grover.ll
	$(CC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) $(SCALE) -D USE_NVGPU -o Grover_driver.o -c Grover_driver.cc
	$(NVCC) -ccbin $(CC) -lnvshmem $(NVCC_FLAGS) $(QIR_BRIDGE_FLAGS) $(NVCC_LIBS) -D USE_NVGPU ../qsharp_wrapper.cu Grover.o Grover_driver.o -o $@

TeleportationSample_nvgpu: TeleportationSample_driver.cc TeleportationSample.ll $(common)
	$(QIRCC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) -o TeleportationSample.o -c TeleportationSample.ll
	$(CC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) $(SCALE) -D USE_NVGPU -o TeleportationSample_driver.o -c TeleportationSample_driver.cc
	$(NVCC) -ccbin $(CC) -lnvshmem $(NVCC_FLAGS) $(QIR_BRIDGE_FLAGS) $(NVCC_LIBS) -D USE_NVGPU ../qsharp_wrapper.cu TeleportationSample.o TeleportationSample_driver.o -o $@

Qrng_nvgpu: Qrng_driver.cc Qrng.ll $(common)
	$(QIRCC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) -o Qrng.o -c Qrng.ll
	$(CC) $(QIRCC_FLAGS) -I$(QIR_BRIDGE_PUBLIC) $(SCALE) -D USE_NVGPU -o Qrng_driver.o -c Qrng_driver.cc
	$(NVCC) -ccbin $(CC) -lnvshmem $(NVCC_FLAGS) $(QIR_BRIDGE_FLAGS) $(NVCC_LIBS) -D USE_NVGPU ../qsharp_wrapper.cu Qrng.o Qrng_driver.o -o $@

clean:
	rm *.o Grover_nvgpu TeleportationSample_nvgpu Qrng_nvgpu

