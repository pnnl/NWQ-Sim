source environment/setup_perlmutter.sh

module load PrgEnv-gnu
module load cmake
module load cpe-cuda
module load cudatoolkit
module unload craype-accel-nvidia80

export CRAYPE_LINK_TYPE=dynamic
export MPICH_GPU_SUPPORT_ENABLED=0
