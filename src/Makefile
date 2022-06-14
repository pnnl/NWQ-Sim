# ---------------------------------------------------------------------------
# NWQsim: Northwest Quantum Circuit Simulation Environment
# ---------------------------------------------------------------------------
# Ang Li, Senior Computer Scientist
# Pacific Northwest National Laboratory(PNNL), U.S.
# Homepage: http://www.angliphd.com
# GitHub repo: http://www.github.com/pnnl/DM-Sim
# PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160
# BSD Lincese.
# ---------------------------------------------------------------------------
# File: Makefile
# Please use CMake for more automatic configuration!!
# ---------------------------------------------------------------------------

include config.mk

CPU_OBJ = test_sv_cpu_omp
#CPU_OBJ = test_sv_cpu_sin
#CPU_OBJ = test_cpu_sin
#CPU_OBJ = test_cpu_omp

#CPU_OBJ = adder_n10_cpu_omp

all: $(CPU_OBJ)

#adder_n10_cpu_sin: adder_n10_cpu_sin.cpp dmsim_cpu_sin.hpp config.h util.h
	#$(CC) $(CC_FLAGS) $(LIBS) adder_n10_cpu_sin.cpp -o $@

test_sv_cpu_sin: test_sv_cpu_sin.cpp svsim_cpu_sin.hpp config.h util.h
	$(CC) $(CC_FLAGS) $(LIBS) test_sv_cpu_sin.cpp -o $@
	#$(CC) $(CC_FLAGS) -mavx512dq -mavx512vl $(LIBS) test_sv_cpu_sin.cpp -o $@


test_sv_cpu_omp: test_sv_cpu_omp.cpp svsim_cpu_omp.hpp nwqsim_cpu_avx512.hpp config.h util.h
	$(CC) $(CC_FLAGS) -mavx512dq -mavx512vl -Xpreprocessor -fopenmp $(LIBS) -lomp test_sv_cpu_omp.cpp -o $@

test_cpu_sin: test_cpu_sin.cpp dmsim_cpu_sin.hpp config.h util.h
	$(CC) $(CC_FLAGS) $(LIBS) test_cpu_sin.cpp -o $@
	#$(CC) $(CC_FLAGS) -mavx512dq -mavx512vl $(LIBS) test_cpu_sin.cpp -o $@

test_cpu_omp: test_cpu_omp.cpp dmsim_cpu_omp.hpp nwqsim_cpu_avx512.hpp config.h util.h
	$(CC) $(CC_FLAGS) -mavx512dq -mavx512vl -Xpreprocessor -fopenmp $(LIBS) -lomp test_cpu_omp.cpp -o $@




adder_n10_cpu_omp: adder_n10_cpu_omp.cpp dmsim_cpu_omp.hpp config.h util.h
	#$(CC) $(CC_FLAGS) -mavx512dq -Xpreprocessor -fopenmp $(LIBS) -lomp adder_n10_cpu_omp.cpp -o $@
	$(CC) $(CC_FLAGS) -mavx512vl -mavx512dq -Xpreprocessor -fopenmp $(LIBS) -lomp adder_n10_cpu_omp.cpp -o $@

clean:
	rm -rf *.o $(OBJ)


