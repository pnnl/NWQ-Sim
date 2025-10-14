# NWQ-VQE

To test ADAPT-VQE on NERSC Perlmutter

## Features

- **High Performance**: GPU-accelerated quantum circuit simulation
- **HPC Scaling**: Leverage multiple CPU/GPU nodes on NERSC Perlmutter

## Building

This module builds together with NWQ-Sim. You need to enable the option of MPI: -DVQE\_ENABLE\_MPI  

```bash
git clone https://github.com/pnnl/NWQ-Sim.git
cd NWQ-Sim
git checkout vqe_memory
git submodule update --init --recursive vqe/nlopt
mkdir build
cd build
source ../environment/setup_perlmutter.sh
module load python
cmake .. -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_CUDA_HOST_COMPILER=CC -DCMAKE_BUILD_TYPE=Release -DVQE_ENABLE_MPI=ON && make -j16
```

## Running

Allocate multi-node GPU jobs on Perlmutter, run on three test cases:


Case-1 (assuming in build dir and using 2 nodes, 1 GPU per node and 128 CPUs per node:
```bash
salloc --nodes 2 --qos interactive -t 120 --constraint gpu --account=m4243
source ../environment/setup_perlmutter.sh
module load python
time srun -n 2 -c 128 --cpu_bind=cores -G 2 --gpu-bind=none vqe/nwq_vqe -b NVGPU -f ../BZDZ-7Orbitals/ducc3_benzene-FrozenCoreCCSD_6Elec_7Orbs.out-xacc -p 6 -v --abstol 1e-6 --maxeval 5000 -o LN_BOBYQA --adapt -ag 1e-3 -am 120
```

Case-2 (assuming in build dir and using 4 nodes, 1 GPU per node and 128 CPUs per node:
```bash
time srun -n 4 -c 128 --cpu_bind=cores -G 4 --gpu-bind=none vqe/nwq_vqe -b NVGPU -f ../BZDZ-9Orbitals/ducc3_benzene-FrozenCoreCCSD_8Elec_9Orbs.out-xacc -p 6 -v --abstol 1e-6 --maxeval 5000 -o LN_BOBYQA --adapt -ag 1e-3 -am 120
```

Case-3 (assuming in build dir and using 4 nodes, 16 MPI tasks, with 1 GPU per task and 32 CPUs per task (suggested for large case, as GPU parts more dominate):
```bash
time srun -N 4 --ntasks-per-node=4 -c 32 -C gpu --gpus-per-task=1 --cpu_bind=cores vqe/nwq_vqe -b NVGPU -f ../H20-11Orbitals/H2O_1.75_Eq_11-Orbitals_DUCC3_H2O-1.75_Eq_DUCC3_10-electrons_11-Orbitals.out-xacc -p 6 -v --abstol 1e-6 --maxeval 5000 -o LN_BOBYQA --adapt -ag 1e-3 -am 120
```

