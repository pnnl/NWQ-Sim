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
# Add the CUDA runtime library

    add_executable(cuda_test cuda_test.cu)
    add_library(CUDA::cudart_static STATIC IMPORTED)
    set_target_properties(CUDA::cudart_static PROPERTIES IMPORTED_LOCATION /opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/11.7/lib64/libcudart_static.a)
    
    set_property(TARGET cuda_test PROPERTY CUDA_ARCHITECTURES 80)
    target_link_libraries(cuda_test PRIVATE CUDA::cudart_static)

    # target_link_libraries(cuda_test PRIVATE CUDA::cudart_static)
endif()
