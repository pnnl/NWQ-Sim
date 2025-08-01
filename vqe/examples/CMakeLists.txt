cmake_minimum_required(VERSION 3.20)

include_directories(BEFORE "include" "../include" ${NLOPT_DIR})
link_directories(BEFORE ${NLOPT_DIR})
link_libraries(vqe)

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
  add_executable(basic_example_cuda "basic_example_cuda.cu")
  set_source_files_properties("basic_example_cuda.cu" PROPERTIES LANGUAGE CUDA)

  target_compile_features(basic_example_cuda PRIVATE cuda_std_17)
  target_compile_definitions(basic_example_cuda PUBLIC CUDA_ENABLED)
  # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        message("NVSMEM FOUND")
        add_executable(basic_example_cuda_mpi "basic_example_cuda_mpi.cu")
        
        set_source_files_properties("basic_example_cuda_mpi.cu" PROPERTIES LANGUAGE CUDA)

        target_compile_features(basic_example_cuda_mpi PRIVATE cuda_std_17)
        target_compile_definitions(basic_example_cuda_mpi PUBLIC CUDA_ENABLED)
        set_target_properties(basic_example_cuda_mpi PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(basic_example_cuda_mpi PUBLIC nvshmem::nvshmem)
        target_compile_definitions(basic_example_cuda_mpi PUBLIC CUDA_MPI_ENABLED)
            target_link_libraries(basic_example_cuda_mpi PUBLIC MPI::MPI_C)
            target_compile_definitions(basic_example_cuda_mpi PUBLIC MPI_ENABLED)
            
            
        # add_executable(cuda_bug_reproduce "cuda_bug_reproduce.cu")
        
        # set_source_files_properties("cuda_bug_reproduce.cu" PROPERTIES LANGUAGE CUDA)

        # target_compile_features(cuda_bug_reproduce PRIVATE cuda_std_17)
        # target_compile_definitions(cuda_bug_reproduce PUBLIC CUDA_ENABLED)
        # set_target_properties(cuda_bug_reproduce PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        # target_link_libraries(cuda_bug_reproduce PUBLIC nvshmem::nvshmem)
        # target_compile_definitions(cuda_bug_reproduce PUBLIC CUDA_MPI_ENABLED)
        #     target_link_libraries(cuda_bug_reproduce PUBLIC MPI::MPI_C)
        #     target_compile_definitions(cuda_bug_reproduce PUBLIC MPI_ENABLED)
    endif()
endif()

# Link against MPI if available
if(MPI_FOUND)
    add_executable(basic_example_mpi "basic_example_mpi.cpp")
    target_link_libraries(basic_example_mpi PUBLIC MPI::MPI_C)
    target_compile_definitions(basic_example_mpi PUBLIC MPI_ENABLED)
endif()

add_executable(basic_example_cpu "basic_example_cpu.cpp")
add_executable(basic_example_adapt_cpu "basic_example_adapt_cpu.cpp")
add_executable(config_example "config_example.cpp")
