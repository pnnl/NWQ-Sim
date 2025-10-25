# ===============================================================================
# VQE Module CMakeLists.txt (new implementation)
# ===============================================================================

# To control if ADAPT VQE use MPI for parallel execution, by default off for daily usage
# should be enabled for HPC.
option(VQE_ENABLE_MPI "Enable MPI for ADAPT_VQE parallelization" OFF)

# Collect source files
set(VQE_SRCS
    src/hamiltonian_parser.cpp
    src/jw_transform.cpp
    src/pauli_term.cpp
    src/environment.cpp
    src/circuit.cpp
    src/ansatz/uccsd.cpp
    src/ansatz/adapt.cpp
    src/backend/gate_builder.cpp
    src/backend/gate_fusion.cpp
    src/backend/statevector_cpu.cpp
    src/execution/adapt_runner.cpp
    src/execution/vqe_runner.cpp
)

set(VQE_GPU_SOURCE src/backend/statevector_gpu.cu)
set(VQE_GPU_STUB   src/backend/statevector_gpu_stub.cpp)

if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    list(APPEND VQE_SRCS ${VQE_GPU_SOURCE})
    set_source_files_properties(${VQE_GPU_SOURCE} PROPERTIES LANGUAGE CUDA)
else()
    list(APPEND VQE_SRCS ${VQE_GPU_STUB})
endif()

# Add nlopt subdirectory first (dependency)
set(_VQE_NLOPT_SOURCE "")
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/nlopt/CMakeLists.txt")
    set(_VQE_NLOPT_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/third_party/nlopt")
elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/nlopt/CMakeLists.txt")
    set(_VQE_NLOPT_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/nlopt")
endif()

if (_VQE_NLOPT_SOURCE)
    if (IS_ABSOLUTE "${_VQE_NLOPT_SOURCE}")
        add_subdirectory(${_VQE_NLOPT_SOURCE} ${CMAKE_CURRENT_BINARY_DIR}/nlopt)
    else()
        add_subdirectory(${_VQE_NLOPT_SOURCE} nlopt)
    endif()
    if (EXISTS "${_VQE_NLOPT_SOURCE}/src/api")
        set(VQE_NLOPT_INCLUDE "${_VQE_NLOPT_SOURCE}/src/api")
    endif()
else()
    message(FATAL_ERROR "NLopt source directory not found. Expected third_party/nlopt or nlopt.")
endif()

if (NOT VQE_NLOPT_INCLUDE)
    set(VQE_NLOPT_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/nlopt/src/api")
endif()

# ===============================================================================
# Define targets
# ===============================================================================

add_library(vqe SHARED ${VQE_SRCS})
add_library(nwqsim "nwqsim_qflow_lib.cpp")
add_executable(nwq_vqe "main.cpp")

# ===============================================================================
# Configure target properties
# ===============================================================================

set(ALL_TARGETS nwqsim vqe nwq_vqe)
foreach(target ${ALL_TARGETS})
    target_compile_features(${target} PRIVATE cxx_std_17)
endforeach()

if(CUDAToolkit_FOUND)
    foreach(target ${ALL_TARGETS})
        target_compile_features(${target} PRIVATE cuda_std_17)
        target_compile_definitions(${target} PUBLIC VQE_ENABLE_CUDA)
    endforeach()
    target_compile_options(vqe PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(vqe PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(vqe PUBLIC NWQSIM_USE_OPENMP)
endif()

if(VQE_ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(vqe PRIVATE MPI::MPI_CXX)
    target_compile_definitions(vqe PRIVATE VQE_ENABLE_MPI)
    target_link_libraries(nwq_vqe PRIVATE MPI::MPI_CXX)
    target_compile_definitions(nwq_vqe PRIVATE VQE_ENABLE_MPI)
else()
    message(STATUS "MPI is OFF for ADAPT_VQE (serial ADAPT loop)")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(vqe PRIVATE
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>:-O3>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>:-ffast-math>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>:-fno-trapping-math>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>:-fno-math-errno>
    )
endif()

# ===============================================================================
# Include directories and link directories
# ===============================================================================

target_include_directories(vqe
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        $<BUILD_INTERFACE:${VQE_NLOPT_INCLUDE}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/nlopt>
)

target_include_directories(nwqsim
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim/vqe>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nwqsim>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${VQE_NLOPT_INCLUDE}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/nlopt>
)

target_include_directories(nwq_vqe PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${VQE_NLOPT_INCLUDE}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

# ===============================================================================
# Basic linking
# ===============================================================================

target_link_libraries(vqe PRIVATE nlopt m)
target_link_libraries(nwqsim PRIVATE vqe nlopt m)
target_link_libraries(nwq_vqe PRIVATE vqe nlopt m)

# ===============================================================================
# Build type specific configurations
# ===============================================================================

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    foreach(target ${ALL_TARGETS})
        target_compile_definitions(${target} PRIVATE NDEBUG)
    endforeach()
endif()

# ===============================================================================
# CUDA support registration (propagate CUDA language to relevant sources)
# ===============================================================================

if(CUDAToolkit_FOUND)
    set(CUDA_SOURCE_FILES
        "src/backend/statevector_gpu.cu"
    )
    foreach(file ${CUDA_SOURCE_FILES})
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
            set_source_files_properties("${file}" PROPERTIES LANGUAGE CUDA)
        endif()
    endforeach()
endif()

# ===============================================================================
# Python bindings (optional)
# ===============================================================================

option(VQE_BUILD_PYTHON "Build Python bindings for the VQE module" ON)

if(VQE_BUILD_PYTHON)
    # Try to find pybind11 via find_package first (for system installations)
    find_package(pybind11 CONFIG QUIET)

    # If not found via find_package, try using add_subdirectory (for git submodule)
    if(NOT pybind11_FOUND)
        set(_VQE_PYBIND11_SOURCE "")
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pybind11/CMakeLists.txt")
            set(_VQE_PYBIND11_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/pybind11")
        elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/pybind11/CMakeLists.txt")
            set(_VQE_PYBIND11_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/third_party/pybind11")
        endif()

        if(_VQE_PYBIND11_SOURCE)
            message(STATUS "pybind11 not found via find_package, using submodule at ${_VQE_PYBIND11_SOURCE}")
            add_subdirectory(${_VQE_PYBIND11_SOURCE} pybind11)
            set(pybind11_FOUND TRUE)
        else()
            message(WARNING "pybind11 not found in find_package or as submodule; skipping VQE Python bindings")
        endif()
    endif()

    if(pybind11_FOUND)
        pybind11_add_module(nwqsim_vqe_core MODULE python/module.cpp)
        target_link_libraries(nwqsim_vqe_core PRIVATE vqe)
        target_include_directories(nwqsim_vqe_core PRIVATE ${VQE_NLOPT_INCLUDE} ${CMAKE_CURRENT_BINARY_DIR}/nlopt)
        set_target_properties(nwqsim_vqe_core PROPERTIES OUTPUT_NAME "_core" PREFIX "")

        # Disable LTO for pybind11 module to work around GCC zstd decompression bug
        # This only affects the Python binding performance, not the core VQE library
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            target_compile_options(nwqsim_vqe_core PRIVATE -fno-lto)
            target_link_options(nwqsim_vqe_core PRIVATE -fno-lto)
        endif()

        # Set proper rpath for the Python module to find libvqe.dylib in the same directory
        if(APPLE)
            set_target_properties(nwqsim_vqe_core PROPERTIES
                INSTALL_RPATH "@loader_path"
                BUILD_WITH_INSTALL_RPATH TRUE)
        else()
            set_target_properties(nwqsim_vqe_core PROPERTIES
                INSTALL_RPATH "$ORIGIN"
                BUILD_WITH_INSTALL_RPATH TRUE)
        endif()

        # Determine installation directory for the Python package
        set(_python_package_subdir "nwqsim_vqe")
        if(DEFINED SKBUILD AND DEFINED SKBUILD_PLATLIB_DIR)
            set(_python_install_dir "${SKBUILD_PLATLIB_DIR}/${_python_package_subdir}")
        else()
            find_package(Python3 COMPONENTS Interpreter REQUIRED)
            execute_process(
                COMMAND "${Python3_EXECUTABLE}" "-c" "import sysconfig;print(sysconfig.get_path('platlib'))"
                OUTPUT_VARIABLE _python_site_package
                OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(NOT _python_site_package)
                message(FATAL_ERROR "Failed to determine Python site-packages directory")
            endif()
            set(_python_install_dir "${_python_site_package}/${_python_package_subdir}")
        endif()

        # Install the Python module and required libraries alongside it
        install(TARGETS nwqsim_vqe_core DESTINATION ${_python_install_dir})
        install(TARGETS vqe DESTINATION ${_python_install_dir})
        install(TARGETS nlopt DESTINATION ${_python_install_dir})

        install(DIRECTORY python/nwqsim_vqe/ DESTINATION ${_python_install_dir} FILES_MATCHING PATTERN "*.py")
    endif()
endif()
