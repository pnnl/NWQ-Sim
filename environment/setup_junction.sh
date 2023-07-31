module load rocm/5.3.0 
module load gcc/10.3.0

export HIP_PATH=/opt/rocm-5.3.0/hip
export PATH=$PATH:$HIP_PATH/bin
export CMAKE_HIP_COMPILER=$HIP_PATH/bin/hipcc