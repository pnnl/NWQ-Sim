# TN SIM 

## Dependencies


### TAMM Dependencies

Baselevel Dependencies
- cmake >= 3.26
- MPI
- C++17 compilter
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

TODO: Add dependencies here


## TAMM Build Instructions

TAMM build instructions can be found here. More direct instructions can be found below specific for Perlmutter.

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

After installing TAMM on Perlmutter you must include the location of the TAMM install directory. 

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


### Personal Computer


## TN Sim Run Instructions

The below commands will allocate a interactive session on Perlmutter across a specified number of nodes. Then run a simulation with TN_TAMM_GPU utilizing all of the GPU's across all of the nodes.

```bash
salloc \
  --nodes=1 \
  --ntasks=5 \
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
  ../environment/perlmutter_bind.sh ./qasm/nwq_qasm -b TN_TAMM_GPU --sim tn --test 3
```



# Recommendatins for Continuation of TAMM NWQ-Sim development

Steps needed for optimization 

## What is Currently Implemented

### Logic for Local and Non-Local 2 qubit gates for MPS 

### Left and Right Environment Orthogonalization

### SVD contraction

## Known errors/bugs/inefficiencies

- SVD is not needed in Left and Right Environment Orthogonalization
- Gates are not executed in parallel
- Measurement is not executed in parallel
- Run time 


## Things that maybe of interest to think about

### 

## Tensor operation parallelization

This also has multiple steps needed for efficient parallelization.

### 1) Levelize sequence of gate operations 

Currently the TAMM NWQ-Sim implementation only schedules a single tensor contraction at a time during simulation. The immediate speed up to this would be to add either in the sim function or a create a new fusion function to schedule all tensor operations that can be performed at once for all sites in the MPS. The best way to do this would be to modify the TAMM code itself. TAMM currently only implents tensor contraction, addition , subtraction etc. SVD is not supported as a operation in the TAMM scheduler. This means that in the current implementation data from TAMM tensors must frequently be taken and inserted into a eigen matrix to perform SVD. 

### 2) TAMM SVD operator implementation

To perform an entire level of operations in a quantum circuit simulation efficiently SVD should be added as a operator that can be added into the scheduler execution graph. Also, current eigen svd is not parallelized across multiple ranks or gpu accelerated. This is a major bottleneck in the current system design. A TAMM svd operator would need to implement this as well.


### 3) Dynamic Tile Size Updates

Currently the tile size is set to a single value. This would also need to become dynamic to the problem size, currently the tile size is set to constant and is not dynamic.



# Final Recommendation

TAMM already implements parallelized tensor contraction distributed in parallel across multiple ranks in a HPC cluster. The tradeoff of rewriting everything from scratch vs updating TAMM code ultimately depends on how much flexibility is needed in the NWQ-Sim TN_Sim design. Implementing the following changes in the TAMM code and the current integration of TAMM into NWQ-Sim would be faster than redesigning everything from scratch


Pro's of TAMM
- MPI, AMD and NVIDIA GPU, memory management already solved in open source PNNL code


Con's of TAMM
- Missing several features for a efficient MPS simulator
- Such as, SVD operator, dynamic tile resizing to deal with Bond dimension
- TAMM seems to have been built in mind to multiply static tensors, a MPS has dynamic tensor sizes depending on the bond dimension
