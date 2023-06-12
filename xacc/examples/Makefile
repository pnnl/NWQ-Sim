include compile_config.mk

XACC_DIR = $(HOME)/.xacc
XACC_LIBS = -I$(XACC_DIR)/include/xacc -I$(XACC_DIR)/include/cppmicroservices4 -I$(XACC_DIR)/include/quantum/gate -L$(XACC_DIR)/lib -lxacc -lxacc-quantum-gate -lCppMicroServices -lxacc-pauli

TARGET = vqe_h2o_adapt.cpp

SRC = $(TARGET) ../nwq_accelerator.cpp 

.PHONY: clean

all: xacc_vqe xacc_qpe xacc_vqe_adapt xacc_qpe_tester

xacc_qpe: qpe.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp 
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(XACC_LIBS) -D USE_CPU qpe.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp  -o $@

xacc_qpe_tester: qpe_tester.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp 
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(XACC_LIBS) -D USE_CPU qpe_tester.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp  -o $@

xacc_vqe: vqe.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp 
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(XACC_LIBS) -D USE_CPU vqe.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp  -o $@

xacc_vqe_adapt: vqe_adapt.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp 
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(XACC_LIBS) -D USE_CPU vqe_adapt.cpp ../backend_runners/sv_omp_runner.cpp ../nwq_accelerator.cpp  -o $@

xacc_adapt_cpu: ../backend_runners/sv_cpu_runner.cpp $(SRC)
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(XACC_LIBS) -D USE_CPU ../backend_runners/sv_cpu_runner.cpp $(SRC) -o $@

xacc_adapt_omp: ../backend_runners/sv_omp_runner.cpp $(SRC)
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(XACC_LIBS) -D USE_CPU -fopenmp ../backend_runners/sv_omp_runner.cpp $(SRC) -o $@

xacc_adapt_nvgpu: ../backend_runners/sv_nvgpu_runner.cu $(SRC)
	$(NVCC_SINGLE) $(NVCC_FLAGS_SINGLE) $(NVCC_LIBS_SINGLE) $(XACC_LIBS) -D USE_NVGPU ../backend_runners/sv_nvgpu_runner.cu $(SRC) -o $@

xacc_adapt_nvgpu_mpi: ../backend_runners/sv_nvgpu_mpi_runner.cu $(SRC)
	$(NVCC_CLUSTER) -ccbin mpicxx -lnvshmem $(NVCC_FLAGS_CLUSTER) $(NVCC_LIBS_CLUSTER) $(SCALE) $(XACC_LIBS) -D USE_NVGPU ../backend_runners/sv_nvgpu_mpi_runner.cu $(SRC) -o $@

clean:
	rm -f xacc_vqe xacc_qpe xacc_vqe_adapt xacc_qpe_tester xacc_adapt_cpu xacc_adapt_omp xacc_adapt_nvgpu xacc_adapt_nvgpu_mpi *.o
