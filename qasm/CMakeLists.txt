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

    set_target_properties(nwq_qasm_nvgpu_mpi PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
                    
    # link libraries
    target_link_libraries(nwq_qasm_nvgpu_mpi PUBLIC nvshmem MPI::MPI_CXX)
endif()
