set(XACC_SOURCE_FILE toy_vqe.cpp nwq_accelerator.cpp)

# Set additional compiler flags for C++
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

add_executable(nwq_xacc ${XACC_SOURCE_FILE})
target_compile_features(nwq_xacc PRIVATE cxx_std_17)

target_link_libraries(nwq_xacc PUBLIC xacc::xacc xacc::quantum_gate)

# Link against OpenMP if available
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(nwq_xacc PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_xacc PUBLIC OMP_ENABLED)
endif()

# Link against MPI if available
if(TARGET MPI::MPI_CXX)
    target_link_libraries(nwq_xacc PUBLIC MPI::MPI_CXX)
    target_compile_definitions(nwq_xacc PUBLIC MPI_ENABLED)
endif()

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties(${XACC_SOURCE_FILE} PROPERTIES LANGUAGE CUDA)

    target_compile_features(nwq_xacc PRIVATE cuda_std_17)
    target_compile_definitions(nwq_xacc PUBLIC CUDA_ENABLED)


    # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        set_target_properties(nwq_xacc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(nwq_xacc PUBLIC nvshmem::nvshmem)
        target_compile_definitions(nwq_xacc PUBLIC CUDA_MPI_ENABLED)
    endif()
endif()

if(HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(${XACC_SOURCE_FILE} PROPERTIES LANGUAGE HIP)

    # target_compile_features(nwq_qasm PRIVATE hip_std_17)
    target_compile_definitions(nwq_xacc PUBLIC HIP_ENABLED)

endif()