# This is for NERSC Perlmutter
make
time srun -C gpu -N 1 -n 1 -c 1 --gpus-per-task=1 --gpu-bind=single:1 ./qasm_dmsim_nvgpu -q ../../../data/openqasm_basis/adder_n10.qasm
