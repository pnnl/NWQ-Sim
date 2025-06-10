# TN SIM 

## Dependencies


### TAMM Dependencies

Baselevel Dependencies
- cmake >= 3.26
- MPI
- C++17 compiler
- CUDA >=11.7
- ROCM >=5.5

Library Dependencies
- GlobalArrays
- HPTT, Librett
- HDF5
- BLAS/LAPACK
- BLAS++ , LAPACK++
- Eigen3, doctest

### iTensor Dependencies

Library Dependencies
- LBLAS, LLAPACK (From LibSci)

## Installation Instructions

General NWQ-Sim installation instructions can be found at:  
[https://github.com/PNNL/NWQ-Sim/blob/main/doc/user_manual.md](https://github.com/PNNL/NWQ-Sim/blob/main/doc/user_manual.md)

This section provides a condensed guide for installing the `tn_sim_tamm_itensor` branch of NWQ-Sim on a local machine or NERSC's Perlmutter cluster.

```bash
git clone https://github.com/PNNL/NWQ-Sim.git -b tn_sim_tamm_itensor
```

## TAMM Build Instructions

TAMM build instructions can be found [here](https://tamm.readthedocs.io/en/latest/install.html). More direct instructions can be found below specific for Perlmutter.

### Perlmutter

Set some temporary install variables.

```bash
export REPO_ROOT_PATH=$HOME/TAMM
export REPO_INSTALL_PATH=$HOME/tamm_install
```

Clone TAMM into Repo path and move to the build directory.

```bash
git clone git@github.com:NWChemEx/TAMM.git $REPO_ROOT_PATH
cd $REPO_ROOT_PATH
mkdir build && cd build
```

Environment setup command, run this when first building and every subsequent time logging on. This will also run the perlmutter setup script.

```bash
source $NWQ_SIM_PATH/environment/setup_tamm_perlmutter.sh
```

Inside the TAMM build directory run the following.

```bash
cmake -DTAMM_ENABLE_CUDA=ON -DGPU_ARCH=80 -DBLIS_CONFIG=generic \
-DCMAKE_INSTALL_PREFIX=$REPO_INSTALL_PATH ..

make -j3
make install
```

After installing TAMM on Perlmutter you must include the location of the TAMM install directory. Modify the following lines in the top level CMakeLists.txt file.

```bash
if ( NOT DEFINED TAMM_DIR )
  set(TAMM_DIR "$ENV{HOME}/PATH/TO/TAMM/INSTALL"
      CACHE PATH "Path to external TAMM installation")
endif()
```

Then the following commands can be run to compile NWQ-Sim for the TN_TAMM_CPU and TN_TAMM_GPU backends.

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DCUDA_ARCH=80
make -j4
```
## iTensor Build Instructions

### Perlmutter

Clone the iTensor repository:

```bash
git clone https://github.com/ITensor/ITensor itensor
```

Edit the `options.mk.sample` file in the `itensor` directory based on your system.  
Ensure BLAS and LAPACK are enabled. For Perlmutter:

```bash
cp ./NWQ-Sim/environment/options.mk ./itensor
source ./NWQ-Sim/environment/setup_perlmutter.sh
cd ./itensor
make
```

```bash
cd ./NWQ-Sim
mkdir build
source ./NWQ-Sim/environment/setup_perlmutter.sh
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCUDA_ARCH=80
make -j4
```

## TN Sim Usage Instructions

### Select Backend

TN-Sim supports three backends. The default is `CPU`.

```bash
./nwq_qasm -q CIRCUIT.QASM --sim tn --backend CPU
./nwq_qasm -q CIRCUIT.QASM --sim tn --backend TN_TAMM_CPU
./nwq_qasm -q CIRCUIT.QASM --sim tn --backend TN_TAMM_GPU
```

### Configure MPS Simulation Parameters

You can configure:

- `--max_dim`: Max number of singular values in SVD (default: 100)
- `--sv_cutoff`: Minimum magnitude threshold for singular values (default: 0.0)

Example:

```bash
./nwq_qasm -q CIRCUIT.QASM --sim tn --max_dim 1000 --sv_cutoff 1e-6
```

### Running TAMM TN-Sim Distributed on Perlmutter

The below commands will allocate a interactive session on Perlmutter across a specified number of nodes. Then run a simulation with TN_TAMM_GPU utilizing all of the GPU's across all of the nodes. The number of tasks to allocate is equal to the number of GPU's per node plus 1 task for the Global Arrays Progress Ranks. The follwing example shows full utilization of 2 nodes on Perlmutter.

```bash
salloc \
  --nodes=2 \
  --ntasks=10 \
  --gpus-per-node=4 \
  --qos=interactive \
  --time=01:00:00 \
  --constraint=gpu \
  --account=mxxxx
```

Inside the build directory run the following to run nwq_qasm with TAMM across multiple nodes and gpus.

```bash
srun -u \
  --cpu_bind=map_cpu:0,16,32,48,64 \
  --mem-bind=map_mem:0,1,2,3,0 \
  --gpus-per-node=4 \
  --ntasks-per-node=5 \
  ../environment/perlmutter_bind.sh ./qasm/nwq_qasm -b TN_TAMM_GPU --sim tn <args>
```

# Recommendations for Continuation of TAMM NWQ-Sim Development

## What is Currently Implemented

- Logic for Local and Non-Local 2-Qubit Gates for MPS
- Full simulation with .qasm files or test
- Measurement function on all sites
- Left and Right Environment Orthogonalization
- Position functoin to set up mixed gauge for left canoncalization on the left and right canonicalization on the right of a specified site
- SVD Contraction utilizing eigen on the host

## Known Errors, Bugs, and Inefficiencies

- SVD is not required in Left and Right Environment Orthogonalization
- Gates are not executed in parallel
- Measurement is not executed in parallel
- Suboptimal run-time performance
- Bond propogation: There is a function in tn_tamm called C2_GATE_NL that is an attempt at this. Currently TAMM is set up to use SWAP gates for non-local 2 qubit gates and is functionally correct.

## Areas of Interest for Further Consideration

### Tensor Operation Parallelization

Parallelization involves multiple stages for efficient execution:

#### 1. Levelize Sequence of Gate Operations

Currently, TAMM NWQ-Sim schedules one tensor contraction at a time. To improve throughput, a new "fusion" function or an enhancement to the existing `sim` function should schedule all tensor operations that can be concurrently executed across all MPS sites. Optimal implementation would require modifying TAMM directly. 

Note: TAMM supports tensor contraction, addition, and subtraction, but not SVD. This forces data transfers between TAMM tensors and Eigen matrices, introducing overhead.

#### 2. TAMM SVD Operator Implementation

To enable efficient execution of entire quantum circuit layers, SVD should be integrated into the TAMM execution graph as a schedulable operator. Additionally, the current Eigen-based SVD implementation is neither parallelized across MPI ranks nor GPU-accelerated. This represents a major bottleneck. A proper TAMM SVD operator should support both.

#### 3. Dynamic Tile Size Updates

The current static tile size configuration is suboptimal. Efficient simulations would require dynamic tile sizing that adapts to changing bond dimensions and tensor shapes.

## Final Recommendation

TAMM already provides a robust foundation for high-performance tensor contraction, including MPI support and GPU memory management for both AMD and NVIDIA platforms. These features are ultimately required to scale a MPS simulation across a HPC cluster such as Perlmutter or Frontier. 

However, TAMM was primarily designed for static tensor contraction and lacks some of the dynamic capabilities required by TN_Sim, such as:

- A native, GPU-accelerated SVD operator
- Support for dynamic tile and tensor sizing, adaptable to changing bond dimensions

To optimize TAMM as a backend for TN_Sim, the following changes are necessary:

1. **Enhance TAMM Core**:
   - Implement a scaled GPU-accelerated SVD operator compatible with the TAMM scheduler
   - Introduce dynamic tile and tensor sizing to eliminate intermediate data copies during simulation

2. **Modify NWQ-Sim Integration**:
   - Adapt the simulation kernel or gate fusion logic to schedule circuit gates using TAMM’s execution model

