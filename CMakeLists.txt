# Specify a minimum CMake version
cmake_minimum_required(VERSION 3.20) 

# Specify the project details
project(MyPackage 
        VERSION 0.0.1 
        DESCRIPTION "NWQSim Qasm" 
        LANGUAGES CXX)

# Set C++ and CUDA standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

#set(CMAKE_C_COMPILER cc)
#set(CMAKE_CXX_COMPILER CC)

# Set the nvshmem home directory
list(APPEND CMAKE_PREFIX_PATH $ENV{HOME}/nvshmem)

# Enable highest optimization for release build
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3")

# Set CUDA architecture
set(CMAKE_CUDA_ARCHITECTURES 80)

# Find required packages
find_package(CUDAToolkit)
find_package(MPI)
find_package(OpenMP)
find_package(NVSHMEM)

# Print package statuses
message(STATUS "CUDA found: ${CUDAToolkit_FOUND}")
message(STATUS "MPI found: ${MPI_FOUND}")
message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
message(STATUS "NVSHMEM found: ${NVSHMEM_FOUND}")

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Subdirectories
add_subdirectory(qasm)
