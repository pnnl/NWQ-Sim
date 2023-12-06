#Environment settings for ORNL Frontier HPC

module use /autofs/nccs-svm1_sw/crusher/amdsw/modules
module load perftools-base/23.05.0
module load cce/16.0.0
module load PrgEnv-cray/8.4.0
module load craype/2.7.21
module load cray-mpich/8.1.26
module load cray-libsci/23.05.1.4
module load craype-accel-amd-gfx90a
module unload cray-pmi/6.1.8
module unload darshan-runtime/3.4.0
module load rocm/5.4.3
module load rocshmem/1.6.1

export cc=gcc
export CC=g++

export LD_LIBRARY_PATH=/opt/rocm-5.4.3/llvm/lib/:$LD_LIBRARY_PATH
export MY_HIP_ARCH=gfx90a
export CRAY_MPICH_ROOTDIR=/opt/cray/pe/mpich/8.1.26
export MPICH_DIR=$CRAY_MPICH_ROOTDIR/ofi/crayclang/14.0
export LDFLAGS="-L${MPICH_DIR}/lib -lmpi ${CRAY_XPMEM_POST_LINK_OPTS} -lxpmem ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a} -I${MPICH_DIR}/include"
export MPICH_GPU_SUPPORT_ENABLED=1
