# This is for NERSC Perlmutter
make
time jsrun -n1 -a1 -g1 -c1 --smpiargs="-gpu"  ./qasm_dmsim_nvgpu -q ../../../data/openqasm_basis/adder_n10.qasm
