#Environment settings for ORNL Frontier HPC
module reset
module load craype-accel-amd-gfx90a
module load cpe/23.12
module load rocm/6.0.0
module load cmake

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
export LDFLAGS="-L${MPICH_DIR}/lib -lmpi ${CRAY_XPMEM_POST_LINK_OPTS} -lxpmem ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a} -I${MPICH_DIR}/include"
export MPICH_GPU_SUPPORT_ENABLED=1
export MY_HIP_ARCH=gfx90a
export rocshmem_DIR=/ccs/home/bpotter/rocshmem
export ROC_SHMEM_MAX_NUM_CONTEXTS=114

export cc=gcc
export CC=g++
