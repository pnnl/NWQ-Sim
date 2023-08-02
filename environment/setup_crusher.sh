#Environment settings for ORNL Crusher HPC
module load rocm
module load cray-mpich
module load craype-accel-amd-gfx90a

export MPICH_GPU_SUPPORT_ENABLED=1
export MY_HIP_ARCH=gfx90a