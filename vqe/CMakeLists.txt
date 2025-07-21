# ===============================================================================
# VQE Module CMakeLists.txt
# ===============================================================================

# Collect source files
set(VQE_SRCS 
    src/jw.cpp              
    src/uccsd.cpp     
    src/utils.cpp
    src/uccsdmin.cpp
    src/singletgsd.cpp  
    src/ansatz_pool.cpp  
    src/hamiltonian.cpp  
)

# Add nlopt subdirectory first (dependency)
add_subdirectory(nlopt)

# ===============================================================================
# Define targets
# ===============================================================================

# VQE base library (internal use)
add_library(vqe SHARED ${VQE_SRCS})

# NWQSim QFlow library object
add_library(nwqsim "nwqsim_qflow_lib.cpp")

# Main VQE executable
add_executable(nwq_vqe "main.cpp")

# QFlow executable 
add_executable(nwq_qflow "qflow.cpp")


# ===============================================================================
# Configure target properties
# ===============================================================================

# Set C++ standard for all targets
set(ALL_TARGETS nwqsim vqe nwq_vqe nwq_qflow)
foreach(target ${ALL_TARGETS})
    target_compile_features(${target} PRIVATE cxx_std_17)
endforeach()

# ===============================================================================
# Include directories and link directories
# ===============================================================================

# Set target-specific include directories
target_include_directories(vqe 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nlopt/src/api>
)

target_include_directories(nwqsim
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim/vqe>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nlopt/src/api>
)

target_include_directories(nwq_vqe PRIVATE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nlopt/src/api>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

target_include_directories(nwq_qflow PRIVATE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/nlopt/src/api>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

# Link directories for internal builds only
target_link_directories(vqe BEFORE PRIVATE "nlopt/build")
target_link_directories(nwqsim BEFORE PRIVATE "nlopt/build")
target_link_directories(nwq_vqe BEFORE PRIVATE "nlopt/build")
target_link_directories(nwq_qflow BEFORE PRIVATE "nlopt/build")

# ===============================================================================
# Basic linking
# ===============================================================================

# Link base VQE library dependencies
target_link_libraries(vqe PRIVATE nlopt m)

# Link NWQSim QFlow library dependencies
target_link_libraries(nwqsim PRIVATE vqe nlopt m)

# Link main VQE executable dependencies
target_link_libraries(nwq_vqe PRIVATE vqe nlopt m)

# Link QFlow executable dependencies
target_link_libraries(nwq_qflow PRIVATE vqe nlopt m)

# ===============================================================================
# Build type specific configurations
# ===============================================================================

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    foreach(target ${ALL_TARGETS})
        target_compile_definitions(${target} PRIVATE NDEBUG)
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
        "nwqsim_qflow_lib.cpp"
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
endif()
