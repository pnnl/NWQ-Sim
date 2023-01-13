# ALCF Theta
make
time aprun -n 4 -N 4 -cc depth -d 1 -j 1 ./qasm_svsim_nvgpu -q ../../../data/openqasm/adder_n10.qasm
