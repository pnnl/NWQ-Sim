# Specify a minimum CMake version
cmake_minimum_required(VERSION 3.20) 


if(APPLE)
    message("Running on Apple platform")
    # Additional configuration for Apple platforms
else()
    message("Not running on Apple platform")
    # Configuration for non-Apple platforms
    set(CMAKE_C_COMPILER "$ENV{cc}")
    set(CMAKE_CXX_COMPILER "$ENV{CC}")
    set(CMAKE_CUDA_HOST_COMPILER mpicxx)

    # Set CUDA architecture
    set(CMAKE_CUDA_ARCHITECTURES "$ENV{MY_CUDA_ARCH}")
    # Set the nvshmem home directory
    list(APPEND CMAKE_PREFIX_PATH $ENV{HOME}/nvshmem)
    set(CMAKE_CUDA_FLAGS_RELEASE "-O3")
endif()


# Specify the project details
project(MyPackage
        VERSION 0.0.1
        DESCRIPTION "NWQSim Qasm"
        LANGUAGES C CXX)

# Set C++ and CUDA standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

if("$ENV{MY_CUDA_ARCH}" STREQUAL "80")
    add_compile_definitions(FP64_TC_AVAILABLE)
endif()


if (NOT XACC_DIR)
  set(XACC_DIR "$ENV{HOME}/.xacc")
endif()

# Find required packages
find_package(CUDAToolkit)
find_package(MPI)
find_package(OpenMP)
find_package(NVSHMEM QUIET)
find_package(XACC QUIET)

if(NVSHMEM_FOUND)
    set(nvshmem_found_str "TRUE")
else()
    set(nvshmem_found_str "FALSE")
endif()

if(XACC_FOUND)
    set(xacc_found_str "TRUE")
else()
    set(xacc_found_str "FALSE")
endif()

# Print package statuses
message(STATUS "CUDA found: ${CUDAToolkit_FOUND}")
message(STATUS "MPI found: ${MPI_FOUND}")
message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
message(STATUS "NVSHMEM found: ${nvshmem_found_str}")
message(STATUS "XACC found: ${xacc_found_str}")

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Subdirectories
add_subdirectory(qasm)

if(XACC_FOUND)
    add_subdirectory(xacc)
endif()
