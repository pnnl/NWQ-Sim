#Environment settings for NERSC Perlmutter HPC
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/qsharp-runtime/src/Qir/Runtime/build/lib/QIR:~/qsharp-runtime/src/Qir/Runtime/build/lib/QSharpCore/:~/qsharp-runtime/src/Qir/Runtime/build/lib/QSharpFoundation/:~/qsharp-runtime/src/Qir/Runtime/build/lib/:~/qsharp-runtime/src/Qir/Runtime/build/lib/Tracer/

module load cmake
module load PrgEnv-nvidia
module load cudatoolkit

export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/comm_libs/mpi/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/nvshmem/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib/:$LD_LIBRARY_PATH

export MY_CUDA_ARCH=80
export NVSHMEM_SYMMETRIC_SIZE=32g

export NVSHMEM_REMOTE_TRANSPORT=libfabric
export NVSHMEM_DISABLE_CUDA_VMM=1
export FI_CXI_DISABLE_HMEM_DEV_REGISTER=1

#XACC https://github.com/eclipse/xacc
export PYTHONPATH=$PYTHONPATH:$HOME/.xacc
export LD_LIBRARY_PATH=$HOME/.xacc/lib:$LD_LIBRARY_PATH

export cc=cc
export CC=CC

# Exports for TAMM

# Link against dynamic libraries
export CRAYPE_LINK_TYPE=dynamic

# Disable MPICH's native GPU support (use OFI instead)
export MPICH_GPU_SUPPORT_ENABLED=0

# Fix one OpenMP thread per MPI rank by default
export OMP_NUM_THREADS=1

# Skip NIC symmetry tests in MPICH‐OFI
export MPICH_OFI_SKIP_NIC_SYMMETRY_TEST=1

# Enable verbose OFI diagnostics (optional; you can unset if too chatty)
export MPICH_OFI_VERBOSE=1
export MPICH_OFI_NIC_VERBOSE=1

# Tell TAMM how many GA progress ranks per node
export GA_NUM_PROGRESS_RANKS_PER_NODE=1
export GA_PROGRESS_RANKS_DISTRIBUTION_PACKED=1

export PPn=4

#Use this when issues with MPI/NVSHMEM
#FI_MR_CACHE_MONITOR=disabled 
#FI_MR_CUDA_CACHE_MONITOR_ENABLED=0
