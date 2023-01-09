# NERSC Perlmutter
make
srun -C gpu -N 1 -n 4 -c 2 --gpus-per-task=1 --gpu-bind=single:1 ./qasm_dmsim_nvgpu -q ../../../data/openqasm_basis/adder_n10.qasm
