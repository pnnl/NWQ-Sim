# Specify a minimum CMake version
cmake_minimum_required(VERSION 3.20) 

# Specify the project details
project(MyPackage 
        VERSION 1.0.0 
        DESCRIPTION "My C++ package" 
        LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable highest optimization for release build
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CUDA_FLAGS_RELEASE "-O3")

# Find CUDA package
find_package(CUDAToolkit)
if(CUDAToolkit_FOUND)
    message(STATUS "CUDA found")
    set(CMAKE_CUDA_ARCHITECTURES 80)
else()
    message(STATUS "CUDA not found")
endif()

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
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
else()
    message(STATUS "OpenMP not found")
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Subdirectories
add_subdirectory(qasm)
