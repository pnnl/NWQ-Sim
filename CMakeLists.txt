cmake_minimum_required(VERSION 3.20)  # Set your version here
project(MyPackage VERSION 1.0.0 DESCRIPTION "My C++ package" LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Specify the directory path to check
set(NVSHMEM_HOME $ENV{HOME}/nvshmem_old)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -O3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# Find MPI
find_package(MPI)

if(MPI_FOUND)
    message(STATUS "MPI found")

    include_directories(${MPI_INCLUDE_PATH})
else()
    message(STATUS "MPI not found")
endif()

# Check if OpenMP is supported
find_package(OpenMP)

if(OpenMP_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "OpenMP not found")
endif()

# Find CUDA package
find_package(CUDA)

if(CUDA_FOUND)
    message(STATUS "CUDA found")
    # Add the CUDA runtime library
    set(CMAKE_CUDA_ARCHITECTURES 80)    
    # Set NVCC flags to enable -O3 optimization
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -arch=sm_80 -m64 -ccbin CC")
else()
    message(STATUS "CUDA not found")
endif()

# Check if NVSHMEM_HOME path exists
if(EXISTS ${NVSHMEM_HOME})
    include_directories(${NVSHMEM_HOME}/include)
    link_directories(${NVSHMEM_HOME}/lib)
    set(LIBS ${LIBS} nvshmem)
    set(NVSHMEM_FOUND TRUE)
    message(STATUS "NVSHMEM found")
else()  # NVSHMEM_HOME path does not exist
  message(WARNING "NVSHMEM_HOME path does not exist. Please set NVSHMEM_HOME to point to the NVSHMEM installation directory.")
endif()


# Set libfabric path
set(LIBFABRIC_PATH "/opt/cray/libfabric/1.15.2.0")

include_directories(${LIBFABRIC_PATH}/include)
link_directories(${LIBFABRIC_PATH}/lib64)
set(LIBS ${LIBS} fabric)

# Include directories
include_directories(include)

# Subdirectories
# add_subdirectory(examples)
add_subdirectory(qasm)
