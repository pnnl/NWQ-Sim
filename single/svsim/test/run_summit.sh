# This is for OLCF Summit
make
time jsrun -n1 -a1 -g1 -c1 --smpiargs="-gpu" ./qasm_svsim_benchmark_nvgpu -a
# ./qasm_svsim_benchmark_nvgpu -a
