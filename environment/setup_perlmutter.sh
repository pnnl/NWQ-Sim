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

#XACC https://github.com/eclipse/xacc
export PYTHONPATH=$PYTHONPATH:$HOME/.xacc
export LD_LIBRARY_PATH=$HOME/.xacc/lib:$LD_LIBRARY_PATH

export cc=cc
export CC=CC

#Use this when issues with MPI/NVSHMEM
#FI_MR_CACHE_MONITOR=disabled 
#FI_MR_CUDA_CACHE_MONITOR_ENABLED=0
