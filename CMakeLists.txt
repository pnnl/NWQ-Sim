# Specify a minimum CMake version
cmake_minimum_required(VERSION 3.20) 

if(DEFINED ENV{cc} AND DEFINED ENV{CC})
    message("Environment variables cc and CC are set")
    set(CMAKE_C_COMPILER "$ENV{cc}")
    set(CMAKE_CXX_COMPILER "$ENV{CC}")
else()
    message("Environment variables cc and/or CC are not set, using default settings")
endif()


# Specify the project details
project(MyPackage
        VERSION 0.0.1
        DESCRIPTION "NWQSim Qasm"
        LANGUAGES C CXX)

# Set local installation prefix to avoid system installation
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/vqe/install" CACHE PATH "Install path prefix" FORCE)

# Set C++ and CUDA standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT XACC_DIR)
  set(XACC_DIR "$ENV{HOME}/.xacc")
endif()

# Find required packages
find_package(CUDAToolkit)
find_package(MPI)
find_package(OpenMP)
find_package(XACC QUIET)
find_package(HIP QUIET)

set(ma_enable "N/A") 
if(CUDAToolkit_FOUND)

    # Set CUDA architecture
    set(CMAKE_CUDA_ARCHITECTURES "$ENV{MY_CUDA_ARCH}")

    # Check for command line option, then environment variable
    if(NOT DEFINED CUDA_ARCH)
        if(DEFINED ENV{MY_CUDA_ARCH})
            set(CMAKE_CUDA_ARCHITECTURES "$ENV{MY_CUDA_ARCH}")
        else()
            message(FATAL_ERROR "CUDA architecture is not set. Please set it manually using cmake -DCUDA_ARCH=N .. or by setting the MY_CUDA_ARCH environment variable.")
        endif()
    else()
        set(CMAKE_CUDA_ARCHITECTURES "${CUDA_ARCH}")
    endif()
    
    # set(CMAKE_CUDA_FLAGS_DEBUG "-g -G")

    set(CMAKE_CUDA_FLAGS_RELEASE "-O3")
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    if("$ENV{MY_CUDA_ARCH}" STREQUAL "80")
        add_compile_definitions(FP64_TC_AVAILABLE)
        set(ma_enable "TRUE")
    else()
        set(ma_enable "FALSE")
    endif()

    # Set the nvshmem home directory
    list(APPEND CMAKE_PREFIX_PATH $ENV{HOME}/nvshmem)

    find_package(NVSHMEM QUIET)

    if(NVSHMEM_FOUND)
        set(nvshmem_found_str "TRUE")

        set(CMAKE_CUDA_HOST_COMPILER mpicxx)
    else()
        set(nvshmem_found_str "FALSE")
    endif()
else()
    set(nvshmem_found_str "N/A")
endif()


if(XACC_FOUND)
    set(xacc_found_str "TRUE")
else()
    set(xacc_found_str "FALSE")
endif()


if(HIP_FOUND)
    set(hip_found_str "TRUE")
    set(CMAKE_HIP_FLAGS_RELEASE "-O3")

    # Check for command line option, then environment variable
    if(NOT DEFINED HIP_ARCH)
        if(DEFINED ENV{MY_HIP_ARCH})
            set(CMAKE_HIP_ARCHITECTURES "$ENV{MY_HIP_ARCH}")
        else()
            message(FATAL_ERROR "HIP architecture is not set. Please set it manually using cmake -DHIP_ARCH=N .. or by setting the MY_HIP_ARCH environment variable.")
        endif()
    else()
        set(CMAKE_HIP_ARCHITECTURES "${HIP_ARCH}")
    endif()
    # set(CMAKE_HIP_ARCHITECTURES "$ENV{MY_HIP_ARCH}")
    # Check for HIP architecture and rocWMMA library
    if("$ENV{MY_HIP_ARCH}" MATCHES "gfx908|gfx90a|gfx940|gfx942")
      # Check if rocWMMA is installed
      find_path(ROCWMMA_INCLUDE_DIR
        NAMES rocwmma/rocwmma.hpp
        PATHS /opt/rocm-5.4.0/include
        NO_DEFAULT_PATH
      )

      if(ROCWMMA_INCLUDE_DIR)
        # Both conditions are met, add the definition
        add_compile_definitions(FP64_HIPMC_AVAILABLE)
        set(ma_enable "TRUE") 
        else()
        # HIP architecture matches, but rocWMMA not found
        set(ma_enable "FALSE, rocWMMA not found ") 
      endif()
    else()
      # HIP architecture doesn't match
      message(WARNING "FP64_HIPMC_AVAILABLE not defined: ")
      set(ma_enable "FALSE, MY_HIP_ARCH does not match required architectures (gfx90a|gfx940|gfx942)") 
    endif()
else()
    set(hip_found_str "FALSE")
endif()

# Print package statuses
message(STATUS "\n==================== NWQSim Environment Search Summary ===================")
message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
message(STATUS "MPI found: ${MPI_FOUND}")
message(STATUS "CUDA found: ${CUDAToolkit_FOUND}")
message(STATUS "HIP found: ${hip_found_str}")
message(STATUS "NVSHMEM found: ${nvshmem_found_str}")
message(STATUS "XACC found: ${xacc_found_str}")
message(STATUS "Matrix Accelerator ENABLED: ${ma_enable}")
message(STATUS "==========================================================================\n")

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Subdirectories
add_subdirectory(qasm)

# Only add VQE subdirectory if nlopt exists - this creates the NWQSimQFlow library
if (EXISTS "${PROJECT_SOURCE_DIR}/vqe/nlopt/CMakeLists.txt")
    add_subdirectory(vqe)
endif()

if(XACC_FOUND)
    add_subdirectory(xacc)
endif()

# Keep the interface library for backward compatibility but don't install the full package
add_library(nwqsim INTERFACE)
add_library(nwqsim::nwqsim ALIAS nwqsim)

target_include_directories(nwqsim INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim>
)
