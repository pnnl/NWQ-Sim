add_executable(nwq_qasm nwq_qasm.cpp)

# Link against OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(nwq_qasm PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_qasm PUBLIC OMP_ENABLED)
endif()

# Link against MPI if available
if(TARGET MPI::MPI_CXX)
    target_link_libraries(nwq_qasm PUBLIC MPI::MPI_CXX)
    target_compile_definitions(nwq_qasm PUBLIC MPI_ENABLED)
endif()

# Only attempt to create CUDA executable if CUDA was found
if(CUDAToolkit_FOUND)
    add_executable(nwq_qasm_nvgpu nwq_qasm_nvgpu.cu)
    target_link_libraries(nwq_qasm_nvgpu CUDA::cudart CUDA::cuda_driver)
endif()

if(CUDAToolkit_FOUND AND TARGET MPI::MPI_CXX)
    add_executable(nwq_qasm_nvgpu_mpi nwq_qasm_nvgpu_mpi.cu)

    target_compile_options(nwq_qasm_nvgpu_mpi PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
        -O3 -arch=sm_80 -m64 -allow-unsupported-compiler -std=c++17 -rdc=true --compiler-options -fPIC -ccbin mpicxx>)

    # include directories
    include_directories("/global/homes/m/mengwang/nvshmem/include")

    # link directories
    link_directories("/global/homes/m/mengwang/nvshmem/lib"
                    "/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64"
                    "/opt/cray/libfabric/1.15.2.0/lib64"
                    "/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64/stubs")

    # link libraries
    target_link_libraries(nwq_qasm_nvgpu_mpi PUBLIC m cuda fabric nvidia nvshmem)
endif()

