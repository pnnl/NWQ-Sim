set(SOURCE_FILE "nwq_qasm.cpp")

add_executable(nwq_qasm nwq_qasm.cpp)
target_compile_features(nwq_qasm PRIVATE cxx_std_17)

# Add linker flag to disable executable stack (fixes CUDA linker warning)
if (NOT APPLE)
    target_link_options(nwq_qasm PRIVATE "-Wl,-z,noexecstack")
endif()

# Link against OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(nwq_qasm PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_qasm PUBLIC OMP_ENABLED)
endif()

if (MPI_FOUND)
  execute_process(
    COMMAND mpicxx -show
    OUTPUT_VARIABLE _mpi_show
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  string(REGEX MATCH "-I[^ ]+" _mpi_include_flag "${_mpi_show}")
  string(REGEX REPLACE "-I" "" _mpi_include_dir "${_mpi_include_flag}")
  message(STATUS "Detected MPI include dir: ${_mpi_include_dir}")

  target_link_libraries(nwq_qasm PUBLIC MPI::MPI_CXX)
  target_compile_definitions(nwq_qasm PUBLIC MPI_ENABLED)
  target_include_directories(nwq_qasm
    PUBLIC
      ${_mpi_include_dir}
      $<$<COMPILE_LANGUAGE:CUDA>:${_mpi_include_dir}>
  )
endif()

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE CUDA)

    target_compile_features(nwq_qasm PRIVATE cuda_std_17)
    # target_link_libraries(nwq_qasm PUBLIC CUDA::cudart)
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

# If iTensor is found, enable itensor TN CPU simulator
if (ITENSOR_FOUND)
    target_link_libraries(nwq_qasm PRIVATE itensor)
    target_compile_definitions(nwq_qasm PUBLIC ITENSOR_ENABLED)
endif()

# If TAMM is found, enable the TN TAMM backend
if (TAMM_FOUND)
    if(CUDAToolkit_FOUND)
        get_target_property(_cudart CUDA::cudart IMPORTED_LOCATION)
        get_filename_component(_cuda_libdir ${_cudart} DIRECTORY)
        list(APPEND CMAKE_BUILD_RPATH ${_cuda_libdir})
    endif()

    if(CUDAToolkit_FOUND)
        target_link_libraries(nwq_qasm PRIVATE
            CUDA::cudart
            CUDA::cublas
            CUDA::cublasLt
            CUDA::cusolver
        )
    endif()

    target_link_libraries(nwq_qasm PRIVATE OpenMP::OpenMP_CXX tamm::tamm)
    target_compile_definitions(nwq_qasm PUBLIC TAMM_ENABLED)
endif()
