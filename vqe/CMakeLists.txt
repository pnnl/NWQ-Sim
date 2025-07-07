# ===============================================================================
# VQE Module CMakeLists.txt
# ===============================================================================

# Collect source files
file(GLOB SOURCEFILES "src/*.cpp")

# Add nlopt subdirectory first (dependency)
add_subdirectory(nlopt)

# ===============================================================================
# Define targets
# ===============================================================================

# VQE base library (internal use)
add_library(vqe SHARED ${SOURCEFILES})

# NWQSim QFlow shared library (for installation)
add_library(nwqsim_qflow_lib SHARED "nwqsim_qflow_lib.cpp")

# Main VQE executable
add_executable(nwq_vqe "main.cpp")

# QFlow executable 
add_executable(nwq_qflow "qflow.cpp")


# ===============================================================================
# Configure target properties
# ===============================================================================

# Set C++ standard for all targets
set(ALL_TARGETS vqe nwqsim_qflow_lib nwq_vqe nwq_qflow)
foreach(target ${ALL_TARGETS})
    target_compile_features(${target} PRIVATE cxx_std_17)
endforeach()

# ===============================================================================
# Include directories and link directories
# ===============================================================================

# Include directories for all targets
include_directories(BEFORE "include" "../include" "nlopt/src/api" "nlopt/build")

# Link directories for internal builds only
target_link_directories(vqe BEFORE PRIVATE "nlopt/build")
target_link_directories(nwqsim_qflow_lib BEFORE PRIVATE "nlopt/build")

# ===============================================================================
# Basic linking
# ===============================================================================

# Link base VQE library dependencies
target_link_libraries(vqe PUBLIC nlopt m)

# Link nwqsim_qflow_lib dependencies
target_link_libraries(nwqsim_qflow_lib PUBLIC vqe nlopt m)

# Link main VQE executable dependencies
target_link_libraries(nwq_vqe PUBLIC vqe)

# Link QFlow executable dependencies
target_link_libraries(nwq_qflow PUBLIC vqe)

# ===============================================================================
# Build type specific configurations
# ===============================================================================

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    foreach(target ${ALL_TARGETS})
        target_compile_definitions(${target} PRIVATE NDEBUG)
    endforeach()
endif()

# ===============================================================================
# Optional dependencies
# ===============================================================================

# OpenMP support
if(OpenMP_FOUND)
    foreach(target ${ALL_TARGETS})
        target_link_libraries(${target} PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(${target} PUBLIC OMP_ENABLED)
    endforeach()
endif()

# MPI support
if(MPI_FOUND)
    foreach(target ${ALL_TARGETS})
        target_link_libraries(${target} PUBLIC MPI::MPI_C)
        target_compile_definitions(${target} PUBLIC MPI_ENABLED)
    endforeach()
endif()

# ===============================================================================
# CUDA support
# ===============================================================================

if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    
    # Set CUDA source file properties
    set(CUDA_SOURCE_FILES 
        "main.cpp" 
        "qflow.cpp" 
        "qflow_py.cpp"
    )
    
    foreach(file ${CUDA_SOURCE_FILES})
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
            set_source_files_properties("${file}" PROPERTIES LANGUAGE CUDA)
        endif()
    endforeach()
    
    # Apply CUDA features and definitions to all targets
    foreach(target ${ALL_TARGETS})
        target_compile_features(${target} PRIVATE cuda_std_17)
        target_compile_definitions(${target} PUBLIC CUDA_ENABLED)
    endforeach()
    
    # NVSHMEM support
    if(NVSHMEM_FOUND AND TARGET nvshmem::nvshmem)
        foreach(target ${ALL_TARGETS})
            set_target_properties(${target} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
            target_link_libraries(${target} PUBLIC nvshmem::nvshmem)
            target_compile_definitions(${target} PUBLIC CUDA_MPI_ENABLED)
        endforeach()
    endif()
endif()

# ===============================================================================
# Add examples subdirectory
# ===============================================================================
# if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/CMakeLists.txt")
#     add_subdirectory("examples")
# endif()

# ===============================================================================
# Installation configuration for nwqsim_qflow_lib
# ===============================================================================

# Set installation directories relative to CMAKE_INSTALL_PREFIX
set(INSTALL_LIB_DIR "lib")
set(INSTALL_INCLUDE_DIR "include")
set(INSTALL_BIN_DIR "bin")

# Install the shared library and its dependency
install(TARGETS nwqsim_qflow_lib vqe
    EXPORT nwqsim_qflow_lib-targets
    RUNTIME DESTINATION ${INSTALL_BIN_DIR}
    LIBRARY DESTINATION ${INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
    INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}
)

# Install the header file (only one needed for external projects)
install(FILES "include/nwqsim_qflow.hpp"
    DESTINATION ${INSTALL_INCLUDE_DIR}
)

# Create and install the CMake package configuration files
install(EXPORT nwqsim_qflow_lib-targets
    FILE nwqsim_qflow_lib-targets.cmake
    NAMESPACE nwqsim::
    DESTINATION ${INSTALL_LIB_DIR}/cmake/nwqsim_qflow_lib
)

# Create a config file for the package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "nwqsim_qflow_lib-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/nwqsim_qflow_lib-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/nwqsim_qflow_lib-config.cmake"
    INSTALL_DESTINATION ${INSTALL_LIB_DIR}/cmake/nwqsim_qflow_lib
)

install(FILES 
    "${CMAKE_CURRENT_BINARY_DIR}/nwqsim_qflow_lib-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/nwqsim_qflow_lib-config-version.cmake"
    DESTINATION ${INSTALL_LIB_DIR}/cmake/nwqsim_qflow_lib
)