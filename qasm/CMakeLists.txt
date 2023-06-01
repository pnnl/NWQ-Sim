# Assume that we have a test file test1.cpp

add_executable(nwq_qasm nwq_qasm.cpp)

# Include directories for test1 target

target_include_directories(nwq_qasm PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(OpenMP_FOUND)
    target_link_libraries(nwq_qasm PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_qasm PUBLIC OMP_ENABLED)
endif()

if(MPI_FOUND)
    target_link_libraries(nwq_qasm PUBLIC MPI::MPI_CXX)
    target_compile_definitions(nwq_qasm PUBLIC MPI_ENABLED)
endif()

if (CUDA_FOUND)
    add_executable(nwq_qasm_nvgpu nwq_qasm_nvgpu.cu)
    # target_link_libraries(cuda_test ${CUDA_LIBRARIES})
endif()
