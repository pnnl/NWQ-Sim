file(GLOB SOURCEFILES "src/*.cpp")

add_executable(nwq_vqe "main.cpp" ${SOURCEFILES})
add_library(vqe SHARED ${SOURCEFILES})

target_compile_features(nwq_vqe PRIVATE cxx_std_17)
target_compile_features(vqe PRIVATE cxx_std_17)

include_directories(BEFORE PRIVATE "include" "../include" "nlopt/src/api" "nlopt/build")
# target_link_directories(vqe BEFORE  PUBLIC "nlopt/build")
target_link_directories(nwq_vqe BEFORE  PUBLIC "nlopt/build")
target_link_libraries(nwq_vqe PUBLIC vqe)
target_link_libraries(vqe PUBLIC nlopt m)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
target_compile_definitions(vqe PUBLIC NDEBUG)
target_compile_definitions(nwq_vqe PUBLIC NDEBUG)
endif()
# Link against OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(nwq_vqe PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_vqe PUBLIC OMP_ENABLED)
    
    target_link_libraries(vqe PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(vqe PUBLIC OMP_ENABLED)
endif()

# Link against MPI if available
if(MPI_FOUND)
    target_link_libraries(nwq_vqe PUBLIC MPI::MPI_CXX)
    target_compile_definitions(nwq_vqe PUBLIC MPI_ENABLED)
    
    target_link_libraries(vqe PUBLIC MPI::MPI_CXX)
    target_compile_definitions(vqe PUBLIC MPI_ENABLED)
endif()

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties(${SOURCEFILES} PROPERTIES LANGUAGE CUDA)

    target_compile_features(nwq_vqe PRIVATE cuda_std_17)
    #target_link_libraries(nwq_vqe PUBLIC CUDA::cudart)
    target_compile_definitions(nwq_vqe PUBLIC CUDA_ENABLED)
    
    target_compile_features(vqe PRIVATE cuda_std_17)
    #target_link_libraries(vqe PUBLIC CUDA::cudart)
    target_compile_definitions(vqe PUBLIC CUDA_ENABLED)

    # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        set_target_properties(nwq_vqe PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(nwq_vqe PUBLIC nvshmem::nvshmem)
        target_compile_definitions(nwq_vqe PUBLIC CUDA_MPI_ENABLED)
        
        set_target_properties(vqe PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(vqe PUBLIC nvshmem::nvshmem)
        target_compile_definitions(vqe PUBLIC CUDA_MPI_ENABLED)
    endif()
endif()

if(HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(${SOURCEFILES} PROPERTIES LANGUAGE HIP)

    # target_compile_features(nwq_vqe PRIVATE hip_std_17)
    target_compile_definitions(nwq_vqe PUBLIC HIP_ENABLED)

    # target_compile_features(vqe PRIVATE hip_std_17)
    target_compile_definitions(vqe PUBLIC HIP_ENABLED)

endif()


add_subdirectory("examples")

