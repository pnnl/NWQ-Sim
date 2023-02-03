# This is for OLCF Crusher
make
time srun -n 1 --ntasks-per-node=1 ./qasm_svsim_benchmark_amdgpu -a
# ./qasm_svsim_benchmark_nvgpu -a
