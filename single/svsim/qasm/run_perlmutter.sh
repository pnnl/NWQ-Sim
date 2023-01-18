# This is for NERSC Perlmutter
# salloc -N 2 -n 8 --qos interactive_ss11 -t 180 --constraint gpu -c 2 -G 8 --gpus-per-task 1 --account=m4142_g
make
time srun -C gpu -N 1 -n 1 -c 1 --gpus-per-task=1 --gpu-bind=single:1 ./qasm_svsim_nvgpu -q ../../../data/openqasm/adder_n10.qasm
