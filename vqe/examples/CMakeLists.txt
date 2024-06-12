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
endif()

# Link against MPI if available
if(MPI_FOUND)
    add_executable(basic_example_mpi "basic_example_mpi.cpp")
    target_link_libraries(basic_example_mpi PUBLIC MPI::MPI_CXX)
    target_compile_definitions(basic_example_mpi PUBLIC MPI_ENABLED)
endif()

add_executable(basic_example_cpu "basic_example_cpu.cpp")
add_executable(config_example "config_example.cpp")
