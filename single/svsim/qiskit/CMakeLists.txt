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

# NVSHMEM requires CMAKE 3.11
MESSAGE(STATUS "CMake version: ${CMAKE_VERSION}")
CMAKE_MINIMUM_REQUIRED(VERSION 3.11...3.18)


SET(PYTHON_EXECUTABLE "/global/homes/a/angli/.conda/envs/py38/bin/python")
SET(PYTHON_INCLUDE_DIRS "/global/homes/a/angli/.conda/envs/py38/include")


OPTION(USE_PYTHON "Build Python Frontend" ON)

# Define backend
OPTION(USE_CPU "Use CPU backend" OFF)
OPTION(USE_NVGPU "Use NVIDIA GPU backend" ON)
OPTION(USE_AMDGPU "Use AMD GPU backend" OFF)

# Define OMP
OPTION(USE_OPENMP "Use OpenMP for scale-up" OFF)

# Define AVX
OPTION(USE_AVX512 "Accelerated by AVX512 on CPU backend" OFF)


project(SVSim VERSION 2.4 LANGUAGES CXX CUDA)

# Set default build type to Release
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: Debug Release"
		FORCE)
ENDIF()

set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "" FORCE)

# Common flags
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC")

SET(CMAKE_C_COMPILER_FORCED TRUE)
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

FIND_PACKAGE(pybind11 CONFIG REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS}/../)

#CONFIGURE_FILE("VERSION.txt" "VERSION.txt" COPYONLY)

IF(USE_CPU)
    ADD_DEFINITIONS(-DUSE_CPU)
    IF(APPLE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -shared -std=c++11 -fPIC -Wl,-undefined,dynamic_lookup ")
    ELSE()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -shared -std=c++11 -fPIC -W ")
    ENDIF()
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

# AVX512
IF(USE_CPU)
    IF(USE_AVX512)
        ADD_DEFINITIONS(-DUSE_AVX512)
        SET(CMAKE_C_FLAGS "-mavx512dq ${CMAKE_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "-mavx512dq ${CMAKE_CXX_FLAGS}")
    ENDIF()
ENDIF()


# NVIDIA GPU Backend
IF(USE_NVGPU)
    ADD_DEFINITIONS(-DUSE_NVGPU)
    FIND_PACKAGE(CUDA REQUIRED)
    SET(CUDA_ARCH_BIN 80)
    IF (NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
        SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_BIN})
    ENDIF()
    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -gencode arch=compute_${CUDA_ARCH_BIN},code=sm_${CUDA_ARCH_BIN} -rdc=true -allow-unsupported-compiler --compiler-options -fPIC")
ENDIF()

IF(USE_PYTHON)
    IF(USE_NVGPU)
        ADD_LIBRARY("svsim" SHARED "svsim_python.cu")
        SET_TARGET_PROPERTIES("svsim" PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
        SET_TARGET_PROPERTIES("svsim" PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        SET_TARGET_PROPERTIES("svsim" PROPERTIES SUFFIX ".so")
        INCLUDE_DIRECTORIES("../src")
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
        INCLUDE_DIRECTORIES(${CUDA_INCLUDES})

        INCLUDE_DIRECTORIES("/global/homes/a/angli/.conda/envs/py38/lib/python3.8/site-packages/mpi4py/include")

        TARGET_LINK_LIBRARIES("svsim" "${CMAKE_EXE_LINKER_FLAGS}")
        TARGET_LINK_LIBRARIES("svsim" ${PYTHON_LIBRARIES})
        TARGET_LINK_LIBRARIES("svsim" ${CUDA_LIBRARIES})
        TARGET_LINK_LIBRARIES("svsim" "-lgcc_s -L/opt/cray/pe/gcc/11.2.0/snos/lib64/")
        TARGET_LINK_LIBRARIES("svsim" "-lgcc -L/opt/cray/pe/gcc/11.2.0/snos/lib64/")
        TARGET_LINK_LIBRARIES("svsim" "-lnvc -lnvcpumath -lnvhpcatm -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/lib")
    ENDIF()
ENDIF()

