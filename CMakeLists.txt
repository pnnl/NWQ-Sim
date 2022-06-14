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
CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

# Define backend
OPTION(USE_CPU "Use CPU backend" ON)
OPTION(USE_NVGPU "Use NVIDIA GPU backend" OFF)
OPTION(USE_AMDGPU "Use AMD GPU backend" OFF)

# Define OMP/MPI
OPTION(USE_OPENMP "Use OpenMP for scale-up" ON)
OPTION(USE_MPI "Use MPI/SHMEM for scale-out" OFF)

# Define AVX
OPTION(USE_AVX512 "Accelerated by AVX512 on CPU backend" OFF)


project(DMSim VERSION 2.4 LANGUAGES CXX)

# Set default build type to Release
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: Debug Release"
		FORCE)
ENDIF()

# Common flags
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")

SET(CMAKE_C_COMPILER_FORCED TRUE)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

FIND_PACKAGE(pybind11 CONFIG REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS}/../)


CONFIGURE_FILE("VERSION.txt" "VERSION.txt" COPYONLY)


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

ADD_LIBRARY("dmsim" SHARED "src/dmsim_python.cpp")
SET_TARGET_PROPERTIES("dmsim" PROPERTIES SUFFIX ".so")
INCLUDE_DIRECTORIES("src")

ADD_LIBRARY("svsim" SHARED "src/svsim_python.cpp")
SET_TARGET_PROPERTIES("svsim" PROPERTIES SUFFIX ".so")
INCLUDE_DIRECTORIES("src")


#IF(STATIC_LINKING)
#IF(NOT APPLE)
#SET(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
#ENDIF()
#SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
#IF(WIN32)
#LIST(INSERT CMAKE_FIND_LIBRARY_SUFFIXES 0 .lib .a)
#ENDIF()
#ENDIF()
#
#IF(APPLE)
#SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "" FORCE)
#ENDIF()
#
#if(STATIC_LINKING)
#SET(Boost_USE_STATIC_LIBS ON)
#ENDIF()
#FIND_PACKAGE(Boost 1.50 COMPONENTS program_options REQUIRED)
#INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
#
#SET(JKU_LIBS
#${MPFR_LIBRARIES}
#${Boost_LIBRARIES})
#
#add_executable(jku_simulator
#src/main.cpp
#src/Simulator.cpp
#src/QASMsimulator.cpp
#src/QASMscanner.cpp
#src/QASMtoken.cpp
#src/QMDDpackage.cpp
#src/textFileUtilities.cpp
#src/qcost.cpp
#src/timing.cpp
#src/QMDDcircuit.cpp
#src/QMDDcomplexD.cpp
#src/QMDDreorder.cpp)
#
#set_target_properties(jku_simulator PROPERTIES
#LINKER_LANGUAGE CXX
#CXX_STANDARD 14)
#
#IF(STATIC_LINKING)
#set_target_properties(jku_simulator PROPERTIES
#LINK_SEARCH_END_STATIC True)
#ENDIF()
#
#
#target_link_libraries(jku_simulator ${JKU_LIBS})
#
#include_directories(src)

