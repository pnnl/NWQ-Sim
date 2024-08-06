file(GLOB SOURCEFILES "src/*.cpp")

add_executable(nwq_vqe "main.cpp")
add_executable(nwq_qflow "qflow.cpp")
# find_path(PYHEADER "Python.h")
# find_package(Python3 COMPONENTS Interpreter Development)
# add_library(nwqflow "qflow_py.cpp")
# set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
# set(PYHEADER ${Python3_EXECUTABLE})
# add_subdirectory(pybind11)
# pybind11_add_module(nwqflow "qflow_py.cpp")
add_subdirectory(nlopt)
add_library(vqe SHARED ${SOURCEFILES})

target_compile_features(nwq_vqe PRIVATE cxx_std_17)
target_compile_features(vqe PRIVATE cxx_std_17)

include_directories(BEFORE "include" "../include" "nlopt/src/api" "nlopt/build")
# target_include_directories(nwqflow BEFORE PRIVATE "../extern/pybind11/include" ${Python3_INCLUDE_DIRS})
target_link_directories(vqe BEFORE  PUBLIC "nlopt/build")
# target_link_directories(nwq_vqe BEFORE  PUBLIC "nlopt/build")
target_link_libraries(nwq_vqe PUBLIC vqe)
target_link_libraries(nwq_qflow PUBLIC vqe)
# target_link_libraries(nwqflow PUBLIC vqe)
target_link_libraries(vqe PUBLIC nlopt m)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
target_compile_definitions(vqe PUBLIC NDEBUG)
target_compile_definitions(nwq_vqe PUBLIC NDEBUG)
endif()
# Link against OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(nwq_vqe PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_vqe PUBLIC OMP_ENABLED)

    target_link_libraries(nwq_qflow PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(nwq_qflow PUBLIC OMP_ENABLED)

    # target_link_libraries(nwqflow PUBLIC OpenMP::OpenMP_CXX)
    # target_compile_definitions(nwqflow PUBLIC OMP_ENABLED)
endif()

# Link against MPI if available
if(MPI_FOUND)
    target_link_libraries(nwq_vqe PUBLIC MPI::MPI_CXX)
    target_compile_definitions(nwq_vqe PUBLIC MPI_ENABLED)
    target_link_libraries(nwq_qflow PUBLIC MPI::MPI_CXX)
    target_compile_definitions(nwq_qflow PUBLIC MPI_ENABLED)

    # target_link_libraries(nwqflow PUBLIC MPI::MPI_CXX)
    # target_compile_definitions(nwqflow PUBLIC MPI_ENABLED)

endif()

# If CUDA Toolkit is available, set the CUDA standard and definitions
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    set_source_files_properties("main.cpp" PROPERTIES LANGUAGE CUDA)
    set_source_files_properties("qflow.cpp" PROPERTIES LANGUAGE CUDA)
    set_source_files_properties("qflow_py.cpp" PROPERTIES LANGUAGE CUDA)

    target_compile_features(nwq_vqe PRIVATE cuda_std_17)
    #target_link_libraries(nwq_vqe PUBLIC CUDA::cudart)
    target_compile_definitions(nwq_vqe PUBLIC CUDA_ENABLED)

    target_compile_features(nwq_qflow PRIVATE cuda_std_17)
    target_compile_definitions(nwq_qflow PUBLIC CUDA_ENABLED)


    # target_compile_features(nwqflow PRIVATE cuda_std_17)
    # target_compile_definitions(nwqflow PUBLIC CUDA_ENABLED)


    # If NVSHMEM is available, link against it and set the definitions
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        set_target_properties(nwq_vqe PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(nwq_vqe PUBLIC nvshmem::nvshmem)
        target_compile_definitions(nwq_vqe PUBLIC CUDA_MPI_ENABLED)

        set_target_properties(nwq_qflow PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(nwq_qflow PUBLIC nvshmem::nvshmem)
        target_compile_definitions(nwq_qflow PUBLIC CUDA_MPI_ENABLED)

        # set_target_properties(nwqflow PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        # target_link_libraries(nwqflow PUBLIC nvshmem::nvshmem)
        # target_compile_definitions(nwqflow PUBLIC CUDA_MPI_ENABLED)
    endif()
endif()


add_subdirectory("examples")

