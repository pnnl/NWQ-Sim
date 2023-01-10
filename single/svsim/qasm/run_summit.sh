# This is for NERSC Perlmutter
make
time jsrun -n1 -a1 -g1 -c1 --smpiargs="-gpu"  ./qasm_svsim_nvgpu -q ../../../data/openqasm/adder_n10.qasm
