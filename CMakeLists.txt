cmake_minimum_required(VERSION 3.15)  # Set your version here
project(MyPackage VERSION 1.0.0 DESCRIPTION "My C++ package" LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(CheckCXXSourceCompiles)

# set(CMAKE_REQUIRED_FLAGS "-mavx512f")
# check_cxx_source_compiles("
#     #include <immintrin.h>
#     int main() {
#         __m512i a = _mm512_set1_epi32(0);
#         return 0;
#     }"
#     HAS_AVX512)
# if(HAS_AVX512)
#     add_definitions(-DUSE_AVX512)
# endif()

# Find MPI
find_package(MPI)

if(MPI_FOUND)
    message(STATUS "MPI found")
else()
    message(STATUS "MPI not found")
endif()

# Check if OpenMP is supported
find_package(OpenMP)

if(OpenMP_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "OpenMP not found")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -m64 -O3 -Wall")


# Include directories
include_directories(include)

# Subdirectories
add_subdirectory(examples)
