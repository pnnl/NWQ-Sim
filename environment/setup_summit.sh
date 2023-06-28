#Environment settings for ORNL Summit HPC
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/work/qsharp-runtime/src/Qir/Runtime/build/lib/QIR:~/work/qsharp-runtime/src/Qir/Runtime/build/lib/QSharpCore/:~/work/qsharp-runtime/src/Qir/Runtime/build/lib/QSharpFoundation/:~/work/qsharp-runtime/src/Qir/Runtime/build/lib/:~/work/qsharp-runtime/src/Qir/Runtime/build/lib/Tracer/:/ccs/home/angli/nvshmem/nvshmem_src_2.6.0-1/build/lib/


export SHMEM_SYMMETRIC_HEAP_SIZE=2G
export NVSHMEM_SYMMETRIC_SIZE=8G
module load cmake
module load gcc/9.3.0
module load cuda/11.4.0
module load gdrcopy/2.3
module load spectrum-mpi

