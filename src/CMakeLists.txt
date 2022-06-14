# ---------------------------------------------------------------------------
# NWQsim: Northwest Quantum Circuit Simulation Environment
# ---------------------------------------------------------------------------
# Ang Li, Senior Computer Scientist
# Pacific Northwest National Laboratory(PNNL), U.S.
# Homepage: http://www.angliphd.com
# GitHub repo: http://www.github.com/pnnl/DM-Sim
# PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160
# BSD Lincese.
# ---------------------------------------------------------------------------
# File: CMakeLists.txt
# ---------------------------------------------------------------------------

# NVSHMEM requires CMAKE 3.11
MESSAGE(STATUS "CMake version: ${CMAKE_VERSION}")
CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

#================================== Configuration =======================================
# Define frontend
OPTION(USE_CXX "Use C++ frontend" OFF)
OPTION(USE_PYTHON "Use Python frontend" ON)
OPTION(USE_QIR "Use QIR (Q#/QDK) frontend" OFF)

# Define backend
OPTION(USE_CPU "Use CPU backend" ON)
OPTION(USE_NVGPU "Use NVIDIA GPU backend" OFF)
OPTION(USE_AMDGPU "Use AMD GPU backend" OFF)

# Define OMP/MPI
OPTION(USE_OPENMP "Use OpenMP for scale-up" ON)
OPTION(USE_MPI "Use MPI/SHMEM for scale-out" OFF)

# AVX512 for CPU backend
OPTION(USE_AVX512 "Accelerated by AVX512 on X86-CPU platforms" OFF)

# ==================================================================
SET(CMAKE_C_COMPILER "/usr/bin/gcc")
SET(CMAKE_CXX_COMPILER "/usr/bin/g++")
SET(CMAKE_C_COMPILER_FORCED TRUE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF (USE_NVGPU)
    PROJECT(DMSim VERSION 2.4 LANGUAGES CXX CUDA)
ELSE()
    PROJECT(DMSim VERSION 2.4 LANGUAGES CXX)
ENDIF()

# Define Clang++ for Q#/QIR
IF (USE_QIR)
    SET(CMAKE_C_COMPILER "/home/lian599/raid/qir/llvm-project/build/bin/clang")
    SET(CMAKE_CXX_COMPILER "/home/lian599/raid/qir/llvm-project/build/bin/clang++")
    SET(QIR_INC_PATH "/home/angli/raid/qir/irina/public_repo/qsharp-runtime/src/QirRuntime/")
    SET(QIR_LIB_PATH "/home/angli/raid/qir/irina/public_repo/qsharp-runtime/src/QirRuntime/build/Linux/Release/")
ENDIF()

# Debug or Release
IF (DEFINED DEBUG)
    MESSAGE("Building DM-Sim in DEBUG mode.")
    ADD_DEFINITIONS(-DDMSIM_DEBUG=${DEBUG})
    SET(CMAKE_BUILD_TYPE Debug)
ELSE()
    MESSAGE("Building DM-Sim in RELEASE mode.")
    SET(CMAKE_BUILD_TYPE Release)
ENDIF()

# Common flags
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")

# Debug or Release
IF(DEFINED DEBUG)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
ELSE()
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ENDIF()

# MPI
IF(USE_MPI)
    ADD_DEFINITIONS(-DUSE_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    INCLUDE_DIRECTORIES(SYSTEM ${MPI_INCLUDE_PATH})
    #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -loshmem")
ENDIF()

# OpenMP
IF(USE_OPENMP)
    ADD_DEFINITIONS(-DUSE_OMP)
    FIND_PACKAGE(OpenMP REQUIRED)
    IF(APPLE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
    ELSE()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    ENDIF()
ENDIF()

# X86 CPU with AVX512
IF(USE_CPU)
    ADD_DEFINITIONS(-DUSE_CPU)
    IF(USE_AVX512)
        ADD_DEFINITIONS(-DUSE_AVX512)
        SET(CMAKE_C_FLAGS "-mavx512f ${CMAKE_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "-mavx512f ${CMAKE_CXX_FLAGS}")
    ENDIF(USE_AVX512)
ENDIF()

# NVIDIA GPU Backend
IF(USE_NVGPU)
    ADD_DEFINITIONS(-DUSE_NVGPU)
    FIND_PACKAGE(CUDA REQUIRED)
    SET(CUDA_ARCH_BIN 70)
    IF (NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
        SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_BIN})
    ENDIF()

    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_${CUDA_ARCH_BIN},code=sm_${CUDA_ARCH_BIN} -rdc=true --compiler-options -fPIC")
    IF(USE_OPENMP)
        SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp")
    ENDIF()
    IF(USE_MPI)
        SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -ccbin mpicxx")
    ENDIF()
ENDIF()

# AMD GPU Backend
IF(USE_AMDGPU)
    ADD_DEFINITIONS(-DUSE_AMDGPU)
    FIND_PACKAGE(hip PATHS /opt/rocm-3.10.0/bin)
    SET(CMAKE_CXX_COMPILER "/opt/rocm-3.10.0/bin/hipcc")
ENDIF()

# Python Frontend
IF(USE_PYTHON)
    FIND_PACKAGE(pybind11 CONFIG REQUIRED)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS}/../)
    IF(USE_NVGPU)
        IF(USE_OPENMP)
            SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-options=' -fopenmp -Wall -shared -std=c++11 -fPIC' ")
        ELSEIF(USE_MPI)
            SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-options=' -Wall -shared -std=c++11 -fPIC' ")
        ELSE()
            SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-options=' -Wall -shared -std=c++11 -fPIC' ")
        ENDIF()
    ENDIF()
    IF(USE_CPU)
        IF(USE_OPENMP)
            IF(APPLE)
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -Wall -shared -std=c++11 -fPIC -Wl,-undefined,dynamic_lookup ")
            ELSE()
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -shared -std=c++11 -fPIC ")
            ENDIF()
        ELSEIF(USE_MPI)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -shared -std=c++11 -fPIC ")
        ELSE()
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -shared -std=c++11 -fPIC -Wl,-undefined,dynamic_lookup ")
        ENDIF()
    ENDIF()
    IF(USE_AMDGPU)
        IF(USE_OPENMP)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -shared -std=c++11 -fPIC ")
        ELSEIF(USE_MPI)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -shared -std=c++11 -fPIC ")
        ELSE()
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -shared -std=c++11 -fPIC ")
        ENDIF()
    ENDIF()
ENDIF()

# Q#/QDK QIR Frontend
IF(USE_QIR)
    INCLUDE_DIRECTORIES("${QIR_INC_PATH}/public/")
    INCLUDE_DIRECTORIES("${QIR_INC_PATH}/test/")
    INCLUDE_DIRECTORIES("src")
    #ADD_SUBDIRECTORY("qir")
ENDIF()


IF(USE_PYTHON)
    IF(USE_NVGPU)
        CONFIGURE_FILE("src/py_wrapper.cpp" "src/py_wrapper.cu" COPYONLY)
        ADD_LIBRARY("dmsim" SHARED "src/py_wrapper.cu")
        SET_TARGET_PROPERTIES("dmsim" PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
        SET_TARGET_PROPERTIES("dmsim" PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        TARGET_LINK_LIBRARIES("dmsim" ${PYTHON_LIBRARIES})
        IF(USE_OPENMP)
            TARGET_LINK_LIBRARIES("dmsim" ${OpenMP_CXX_FLAGS})
        ELSEIF(USE_MPI)
            TARGET_LINK_LIBRARIES("dmsim" ${MPI_LIBRARIES})
        ENDIF()
    ENDIF()
    IF(USE_CPU)
        ADD_LIBRARY("dmsim" SHARED "python_interface.cpp")
        IF(USE_OPENMP)
            IF(APPLE)
                TARGET_LINK_LIBRARIES("dmsim" "-lomp")
            ENDIF()
        ELSEIF(USE_MPI)
            TARGET_LINK_LIBRARIES("dmsim" ${MPI_LIBRARIES})
        ENDIF()
    ENDIF()
    IF(USE_AMDGPU)
        ADD_LIBRARY("dmsim" SHARED "src/py_wrapper.cpp")
        IF(USE_MPI)
            TARGET_LINK_LIBRARIES("dmsim" ${MPI_LIBRARIES})
        ENDIF()
    ENDIF()
ENDIF()

INCLUDE_DIRECTORIES("src")

IF(USE_CXX)
    ADD_SUBDIRECTORY("example")
    ADD_SUBDIRECTORY("unit_test")
ENDIF()

#ADD_SUBDIRECTORY("qasmbench")
#ADD_SUBDIRECTORY("qir")

