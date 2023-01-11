# OLCF Summit
make
jsrun -n2 -a1 -g1 -c1 --smpiargs="-gpu"  ./qasm_svsim_benchmark_nvgpu -a
