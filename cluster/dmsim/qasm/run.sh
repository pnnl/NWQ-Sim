# NERSC Perlmutter
make
srun -C gpu -N 2 -n 8 -c 2 --gpus-per-task=1 --gpu-bind=single:1 ./qasm_svsim_nvgpu -q ../../../data/openqasm/adder_n10.qasm
