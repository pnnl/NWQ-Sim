# Specify a minimum CMake version
cmake_minimum_required(VERSION 3.20) 

if(DEFINED ENV{cc} AND DEFINED ENV{CC})
    message("Environment variables cc and CC are set")
    set(CMAKE_C_COMPILER "$ENV{cc}")
    set(CMAKE_CXX_COMPILER "$ENV{CC}")
else()
    message("Environment variables cc and/or CC are not set, using default settings")
endif()

option(NWQSIM_ENABLE_VQE "Build VQE components (requires NLopt)" OFF)


# Specify the project details
project(NWQ-Sim
    VERSION 1.0.0
    DESCRIPTION "NWQ-Sim Quantum Simulation Toolkit"
    LANGUAGES C CXX)

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(nwqsim-utils)

# Set C++ and CUDA standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
nwqsim_option(CMAKE_BUILD_TYPE Release)

# Find required packages
find_package(CUDAToolkit)
find_package(MPI)
find_package(OpenMP)
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
message(STATUS "Matrix Accelerator ENABLED: ${ma_enable}")
message(STATUS "Build VQE component: ${NWQSIM_ENABLE_VQE}")
message(STATUS "==========================================================================\n")


# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Subdirectories
add_subdirectory(qasm)

if(NWQSIM_ENABLE_VQE)
    message(STATUS "NWQSIM_ENABLE_VQE=ON (will look for NLopt)")
    message(STATUS "If NLopt is not found, use 'git submodule update --init --recursive vqe/nlopt'")
    add_subdirectory(vqe)
endif()

set(NWQSIM_INSTALL_TARGETS nwqsim)
# Create nwqsim library - if VQE exists, interface otherwise
if (TARGET nwqsim)
    target_include_directories(nwqsim PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim/vqe>
    )
    list(APPEND NWQSIM_INSTALL_TARGETS vqe nlopt)
    install(FILES
        "vqe/include/nwqsim_qflow.hpp"
        "vqe/include/vqe_options.hpp"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nwqsim/vqe)
    install(FILES
        "vqe/include/execution/adapt_runner.hpp"
        "vqe/include/execution/vqe_runner.hpp"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nwqsim/vqe/execution)
else()
    add_library(nwqsim INTERFACE)
    target_include_directories(nwqsim INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim>
    )
endif()

add_library(nwqsim::nwqsim ALIAS nwqsim)


# Install
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/nwqsim)

# Install main nwqsim library

install(TARGETS ${NWQSIM_INSTALL_TARGETS}
    EXPORT nwqsim-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set(NWQSIM_HEADERS
    include/backendManager.hpp
    include/circuit.hpp
    include/config.hpp
    include/gate.hpp
    include/nwq_util.hpp
    include/state.hpp
)

install(FILES ${NWQSIM_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nwqsim)
install(DIRECTORY include/svsim DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nwqsim)
install(DIRECTORY include/dmsim DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nwqsim)
install(DIRECTORY include/circuit_pass DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nwqsim)

# VQE headers are handled by the vqe/CMakeLists.txt installation

#Export the targets to a script
install(EXPORT nwqsim-targets
    FILE
        nwqsim-targets.cmake
    NAMESPACE
        nwqsim::
    DESTINATION
        ${INSTALL_CONFIGDIR}
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/nwqsim-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/nwqsim-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/nwqsim-config.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
)

#Install the config, configversion and custom find modules
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nwqsim-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/nwqsim-config-version.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

#Export build tree
export(EXPORT nwqsim-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/nwqsim-targets.cmake
    NAMESPACE nwqsim::)

#Register package in user's package registry
# export(PACKAGE nwqsim)
