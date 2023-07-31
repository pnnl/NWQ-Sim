###################### AMD GPU Configuration #######################
HIPCC = /opt/rocm-5.3.0/bin/hipcc
HIPCC_FLAGS = -O3 -std=c++17 
HIPCC_LIBS = -lm 
#######################################################################


qasm_svsim_amdgpu: nwq_qasm_amd.cpp 
	$(HIPCC) $(HIPCC_FLAGS) $(HIPCC_LIBS) -I../include nwq_qasm_amd.cpp -o $@