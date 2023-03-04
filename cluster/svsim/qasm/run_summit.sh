# OLCF Summit
make
jsrun -n8 -a1 -g1 -c2 -r4 --smpiargs="-gpu"  ./qasm_svsim_nvgpu -q ../../../data/openqasm/adder_n10.qasm
