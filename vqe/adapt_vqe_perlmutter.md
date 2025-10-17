# NWQ-VQE (NERSC Perlmutter)

To test ADAPT-VQE on NERSC [Perlmutter](https://docs.nersc.gov/systems/perlmutter/architecture/). Will have a similar one for OLCF Frontier.

## Features

- **High Performance**: GPU-accelerated quantum circuit simulation
- **HPC Scaling**: Leverage multiple CPU/GPU nodes on NERSC Perlmutter

## Building

This module builds together with NWQ-Sim. You need to enable the option of MPI: -DVQE\_ENABLE\_MPI by using the following cmake command or adjust the option of option(VQE_ENABLE_MPI "Enable MPI for ADAPT_VQE parallelization" OFF) from OFF to ON in the cmake file "NWQ-Sim/vqe/CMakeLists.txt". 

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


### Case-1 (assuming in build dir and using 2 nodes, 1 GPU per node and 128 CPUs per node:
```bash
salloc --nodes 2 --qos interactive -t 120 --constraint gpu --account=m4243
source ../environment/setup_perlmutter.sh
module load python
time srun -n 2 -c 128 --cpu_bind=cores -G 2 --gpu-bind=none vqe/nwq_vqe -b NVGPU -f ../BZDZ-7Orbitals/ducc3_benzene-FrozenCoreCCSD_6Elec_7Orbs.out-xacc -p 6 -v --abstol 1e-6 --maxeval 5000 -o LN_BOBYQA --adapt -ag 1e-3 -am 120
```

### Case-2 (assuming in build dir and using 4 nodes, 1 GPU per node and 128 CPUs per node:
```bash
salloc --nodes 4 --qos interactive -t 120 --constraint gpu --account=m4243
source ../environment/setup_perlmutter.sh
module load python
time srun -n 4 -c 128 --cpu_bind=cores -G 4 --gpu-bind=none vqe/nwq_vqe -b NVGPU -f ../BZDZ-9Orbitals/ducc3_benzene-FrozenCoreCCSD_8Elec_9Orbs.out-xacc -p 6 -v --abstol 1e-6 --maxeval 5000 -o LN_BOBYQA --adapt -ag 1e-3 -am 120
```

### Case-3 (assuming in build dir and using 4 nodes, 16 MPI tasks, with 1 GPU per task and 32 CPUs per task (suggested for large case, as GPU parts more dominate, try me!):
```bash
salloc --nodes 4 --qos interactive -t 120 --constraint gpu --account=m4243
source ../environment/setup_perlmutter.sh
module load python
time srun -N 4 --ntasks-per-node=4 -c 32 -C gpu --gpus-per-task=1 --cpu_bind=cores vqe/nwq_vqe -b NVGPU -f ../H20-11Orbitals/H2O_1.75_Eq_11-Orbitals_DUCC3_H2O-1.75_Eq_DUCC3_10-electrons_11-Orbitals.out-xacc -p 6 -v --abstol 1e-6 --maxeval 5000 -o LN_BOBYQA --adapt -ag 1e-3 -am 120
```

## Examplar job script for Perlmutter
```bash
#!/bin/bash

#SBATCH -N 16
#SBATCH -G 64
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -J h20_11a
#SBATCH --mail-user=xxx.xx@pnnl.gov
#SBATCH -o printouts/out_%x_%j.txt
#SBATCH -e printouts/err_%x_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --time-min=32:00:00
#SBATCH -t 10:00:00
#SBATCH -A m4243


export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# VQE parameters
NWQSIM_FOLDER="NWQ-Sim"
BACKEND="NVGPU"

HAMILTONIAN_PATH=H2O_1.75_Eq_11-Orbitals_DUCC3_H2O-1.75_Eq_DUCC3_10-electrons_11-Orbitals.out-xacc
NUM_PARTICLES=10

OPTIMIZER="LN_BOBYQA" 
MAX_ITERATIONS=30000
MAX_ADAITERS=100 
TOLERANCE=1e-6
UB=0.9
LB=-0.9
ADAPT_GRADIENT_THRESHOLD=1e-3
MAX_SECS=3600

echo "Starting ADAPT-VQE simulation on Perlmutter GPU node (GPU backend)"
echo "Job ID: ${SLURM_JOBID}"
echo "Node: ${SLURM_NODELIST}"
echo "GPUs: ${SLURM_GPUS}"
echo "Tasks: ${SLURM_NTASKS}"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "Backend: ${BACKEND}"
echo "Timestamp: $(date)"
echo "================================================"

source ~/${NWQSIM_FOLDER}/environment/setup_perlmutter.sh
module load python

time srun -n 64 -c 32 --cpu_bind=cores -G 64 --gpus-per-task=1 --cpu_bind=cores ${NWQSIM_FOLDER}/build/vqe/nwq_vqe \
    -b ${BACKEND} \
    -f ${HAMILTONIAN_PATH} \
    -p ${NUM_PARTICLES} \
    -v \
    --abstol ${TOLERANCE} \
    -ub ${UB} \
    -lb ${LB} \
    --maxeval ${MAX_ITERATIONS} \
    -o ${OPTIMIZER} \
    --sym 3 \
    --adapt \
    -ag ${ADAPT_GRADIENT_THRESHOLD} \
    -am ${MAX_ADAITERS} \
    --as



echo "================================================"
echo "VQE simulation completed"
echo "End timestamp: $(date)"
```

