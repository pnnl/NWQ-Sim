set(SOURCE_FILE "nwq_qasm.cpp")
find_package(TAMM REQUIRED)

# Tell CMake to use CC for linking C++ executables:
set(CMAKE_CXX_LINK_EXECUTABLE
    "${CMAKE_CXX_COMPILER} <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>"
)

# ── Main executable: nwq_qasm ───────────────────────────────────────────────
add_executable(nwq_qasm ${SOURCE_FILE})
target_compile_features(nwq_qasm PRIVATE cxx_std_17)

# OpenMP
if(OpenMP_FOUND)
  target_link_libraries(nwq_qasm PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(nwq_qasm PUBLIC OMP_ENABLED)
endif()

# MPI
if(MPI_FOUND)
  target_link_libraries(nwq_qasm PUBLIC MPI::MPI_CXX)
  target_compile_definitions(nwq_qasm PUBLIC MPI_ENABLED)
endif()

# CUDA
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE CUDA)

    target_compile_features(nwq_qasm PRIVATE cuda_std_17)
    target_compile_definitions(nwq_qasm PUBLIC CUDA_ENABLED)

    # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        set_target_properties(nwq_qasm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(nwq_qasm PUBLIC nvshmem::nvshmem)
        target_compile_definitions(nwq_qasm PUBLIC CUDA_MPI_ENABLED)
    endif()
endif()

# If HIP is available, enable HIP support
if(HIP_FOUND)
  enable_language(HIP)
  set_source_files_properties(${SOURCE_FILE} PROPERTIES LANGUAGE HIP)
  target_compile_definitions(nwq_qasm PUBLIC HIP_ENABLED)
endif()

# —————————————————————————————————————————————————————
# Link ITensor and Cray LibSci BLAS/LAPACK
target_link_libraries(nwq_qasm
  PRIVATE
    itensor         # imported ITensor library
    sci_nvidia_mpi_mp
    sci_nvidia_mp
)

target_link_libraries(nwq_qasm PUBLIC tamm::tamm)

# ── Smoke‐test executable: tamm_test ────────────────────────────────────────
add_executable(tamm_test tamm_test.cpp)
target_compile_features(tamm_test PRIVATE cxx_std_17)

target_link_libraries(tamm_test PUBLIC
  MPI::MPI_CXX
  OpenMP::OpenMP_CXX
  tamm::tamm
)
