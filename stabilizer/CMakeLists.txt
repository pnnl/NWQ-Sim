set(SOURCE_FILE "circuit_test.cpp")

add_executable(circuit_test circuit_test.cpp)
target_compile_features(circuit_test PRIVATE cxx_std_17)



# Link against OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(circuit_test PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(circuit_test PUBLIC OMP_ENABLED)
endif()

# Link against MPI if available
if(MPI_FOUND)
    target_link_libraries(circuit_test PUBLIC MPI::MPI_CXX)
    target_compile_definitions(circuit_test PUBLIC MPI_ENABLED)
endif()

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE CUDA)

    # Define CUDA_ENABLED for the C++ executable
    # target_compile_definitions(backendManager PRIVATE CUDA_ENABLED)
    set(ENV{MY_CUDA_ARCH} "80")

    target_compile_features(circuit_test PRIVATE cuda_std_17)
    #target_link_libraries(circuit_test PUBLIC CUDA::cudart)
    target_compile_definitions(circuit_test PUBLIC CUDA_ENABLED)

    # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        set_target_properties(circuit_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(circuit_test PUBLIC nvshmem::nvshmem)
        target_compile_definitions(circuit_test PUBLIC CUDA_MPI_ENABLED)
    endif()
endif()

if(HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE HIP)

    # target_compile_features(circuit_test PRIVATE hip_std_17)
    target_compile_definitions(circuit_test PUBLIC HIP_ENABLED)

endif()




set(SOURCE_FILE "surface_code.cpp")

add_executable(surface_code surface_code.cpp)

target_compile_features(surface_code PRIVATE cxx_std_17)



# Link against OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(surface_code PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(surface_code PUBLIC OMP_ENABLED)
endif()

# Link against MPI if available
if(MPI_FOUND)
    target_link_libraries(surface_code PUBLIC MPI::MPI_CXX)
    target_compile_definitions(surface_code PUBLIC MPI_ENABLED)
endif()

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE CUDA)

    # Define CUDA_ENABLED for the C++ executable
    # target_compile_definitions(backendManager PRIVATE CUDA_ENABLED)
    set(ENV{MY_CUDA_ARCH} "80")

    target_compile_features(surface_code PRIVATE cuda_std_17)
    #target_link_libraries(surface_code PUBLIC CUDA::cudart)
    target_compile_definitions(surface_code PUBLIC CUDA_ENABLED)

    # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        set_target_properties(surface_code PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(surface_code PUBLIC nvshmem::nvshmem)
        target_compile_definitions(surface_code PUBLIC CUDA_MPI_ENABLED)
    endif()
endif()

if(HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE HIP)

    # target_compile_features(surface_code PRIVATE hip_std_17)
    target_compile_definitions(surface_code PUBLIC HIP_ENABLED)

endif()


if(HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE HIP)

    # target_compile_features(damping_channel PRIVATE hip_std_17)
    target_compile_definitions(damping_channel PUBLIC HIP_ENABLED)

endif()