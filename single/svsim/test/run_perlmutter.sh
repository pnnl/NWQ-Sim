# This is for NERSC Perlmutter
time srun -C gpu -N 1 -n 1 -c 1 --gpus-per-task=1 --gpu-bind=single:1 ./qasm_svsim_benchmark_nvgpu -a
# ./qasm_svsim_benchmark_nvgpu -a
