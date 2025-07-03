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

# Main VQE executable
add_executable(nwq_vqe "main.cpp")

# QFlow executable 
add_executable(nwq_qflow "qflow.cpp")

# ExaChem QFlow shared library (for external use)
add_library(exachem_qflow_lib SHARED "exachem_qflow_lib.cpp" ${SOURCEFILES})

# Configure build interface include directories only
target_include_directories(exachem_qflow_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

# ===============================================================================
# Configure target properties
# ===============================================================================

# Set C++ standard for all targets
set(ALL_TARGETS vqe nwq_vqe nwq_qflow exachem_qflow_lib)
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
target_link_directories(exachem_qflow_lib BEFORE PRIVATE "nlopt/build")

# ===============================================================================
# Basic linking
# ===============================================================================

# Link base VQE library dependencies
target_link_libraries(vqe PUBLIC nlopt m)

# Create the ExaChem QFlow library with proper dependencies
target_link_libraries(exachem_qflow_lib PUBLIC nlopt m)
# Note: We include VQE source files directly rather than linking to avoid export issues

# Link executables
target_link_libraries(nwq_vqe PUBLIC vqe)
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
        "exachem_qflow_lib.cpp"
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
# Installation configuration for ExaChem QFlow library
# ===============================================================================

# Set installation directory within the project - override CMAKE_INSTALL_PREFIX
set(EXACHEM_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/install")
set(CMAKE_INSTALL_PREFIX "${EXACHEM_INSTALL_DIR}" CACHE PATH "Install path prefix" FORCE)

# Install the ExaChem QFlow shared library
install(TARGETS exachem_qflow_lib
    EXPORT ExaChemQFlowTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install the header file
install(FILES "include/exachem_qflow.hpp" 
    DESTINATION include
)

# Set the install interface include directories after installation
set_target_properties(exachem_qflow_lib PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$<INSTALL_INTERFACE:include>"
)

# Install export targets
install(EXPORT ExaChemQFlowTargets
    FILE ExaChemQFlowTargets.cmake
    NAMESPACE ExaChemQFlow::
    DESTINATION cmake
)

# Create and install CMake config file for external packages
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/ExaChemQFlowConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ExaChemQFlowConfig.cmake"
    @ONLY
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ExaChemQFlowConfig.cmake"
    DESTINATION cmake
)

# ===============================================================================
# Add examples subdirectory
# ===============================================================================

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/CMakeLists.txt")
    add_subdirectory("examples")
endif()