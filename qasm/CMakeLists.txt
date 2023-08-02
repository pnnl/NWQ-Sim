set(SOURCE_FILE "nwq_qasm.cpp")

add_executable(nwq_qasm nwq_qasm.cpp)
target_compile_features(nwq_qasm PRIVATE cxx_std_17)

# Link against OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(nwq_qasm PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_qasm PUBLIC OMP_ENABLED)
endif()

# Link against MPI if available
if(MPI_FOUND)
    target_link_libraries(nwq_qasm PUBLIC MPI::MPI_CXX)
    target_compile_definitions(nwq_qasm PUBLIC MPI_ENABLED)
endif()

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE CUDA)

    target_compile_features(nwq_qasm PRIVATE cuda_std_17)
    #target_link_libraries(nwq_qasm PUBLIC CUDA::cudart)
    target_compile_definitions(nwq_qasm PUBLIC CUDA_ENABLED)

    # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        set_target_properties(nwq_qasm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(nwq_qasm PUBLIC nvshmem::nvshmem)
        target_compile_definitions(nwq_qasm PUBLIC CUDA_MPI_ENABLED)
    endif()
endif()

if(HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE HIP)

    # target_compile_features(nwq_qasm PRIVATE hip_std_17)
    target_compile_definitions(nwq_qasm PUBLIC HIP_ENABLED)

endif()