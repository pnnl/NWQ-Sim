#Environment settings for NERSC Perlmutter HPC
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/qsharp-runtime/src/Qir/Runtime/build/lib/QIR:~/qsharp-runtime/src/Qir/Runtime/build/lib/QSharpCore/:~/qsharp-runtime/src/Qir/Runtime/build/lib/QSharpFoundation/:~/qsharp-runtime/src/Qir/Runtime/build/lib/:~/qsharp-runtime/src/Qir/Runtime/build/lib/Tracer/

module load cmake
module load PrgEnv-nvidia
module load cudatoolkit

export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/comm_libs/mpi/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/nvshmem/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib/:$LD_LIBRARY_PATH



#export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/comm_libs/nvshmem/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/22.5/compilers/lib/:$LD_LIBRARY_PATH

export NVSHMEM_SYMMETRIC_SIZE=32g

export NVSHMEM_REMOTE_TRANSPORT=libfabric
export NVSHMEM_DISABLE_CUDA_VMM=1

#XACC https://github.com/eclipse/xacc
export PYTHONPATH=$PYTHONPATH:$HOME/.xacc
export LD_LIBRARY_PATH=$HOME/.xacc/lib:$LD_LIBRARY_PATH

#module unload cgpu
#module load gcc/9.3.0
#module load cudatoolkit/21.3_11.2
#module load nersc-easybuild/21.12 
#module load GDRCopy/2.1-GCCcore-10.2.0-CUDA-11.1.1
#module load OpenMPI/4.0.5-gcccuda-2020b