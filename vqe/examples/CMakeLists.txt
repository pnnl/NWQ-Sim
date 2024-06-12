cmake_minimum_required(VERSION 3.20)

include_directories(BEFORE PRIVATE "include" "../include" ${NLOPT_DIR})
link_directories(BEFORE  PRIVATE ${NLOPT_DIR})
link_libraries(vqe)

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    # Set CUDA architecture
    set(CMAKE_CUDA_ARCHITECTURES "$ENV{MY_CUDA_ARCH}")

    # Check for command line option, then environment variable
    if(NOT DEFINED CUDA_ARCH)
        if(DEFINED ENV{MY_CUDA_ARCH})
            set(CMAKE_CUDA_ARCHITECTURES "$ENV{MY_CUDA_ARCH}")
        else()
            message(FATAL_ERROR "CUDA architecture is not set. Please set it manually using cmake -DCUDA_ARCH=N .. or by setting the MY_CUDA_ARCH environment variable.")
        endif()
    else()
        set(CMAKE_CUDA_ARCHITECTURES "${CUDA_ARCH}")
    endif()
  add_executable(basic_example_cuda "basic_example_cuda.cu" ${SOURCE})

  enable_language(CUDA)
  set_source_files_properties("basic_example_cuda.cu" PROPERTIES LANGUAGE CUDA)

  target_compile_features(basic_example_cuda PRIVATE cuda_std_17)
  target_compile_definitions(basic_example_cuda PUBLIC CUDA_ENABLED)

  # If NVSHMEM is available, link against it and set the definitions
  if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
      set_target_properties(basic_example_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
      target_link_libraries(basic_example_cuda PUBLIC nvshmem::nvshmem)
      target_compile_definitions(basic_example_cuda PUBLIC CUDA_MPI_ENABLED)
  endif()
endif()

# Link against MPI if available
if(MPI_FOUND)
    add_executable(basic_example_mpi "basic_example_mpi.cpp" ${SOURCE})
    target_link_libraries(basic_example_mpi PUBLIC MPI::MPI_CXX)
    target_compile_definitions(basic_example_mpi PUBLIC MPI_ENABLED)
endif()

add_executable(basic_example_cpu "basic_example_cpu.cpp" ${SOURCE})
add_executable(config_example "config_example.cpp" ${SOURCE})
