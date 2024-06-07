cmake_minimum_required(VERSION 3.20)

if(DEFINED ENV{cc} AND DEFINED ENV{CC})
    message("Environment variables cc and CC are set")
    set(CMAKE_C_COMPILER "$ENV{cc}")
    set(CMAKE_CXX_COMPILER "$ENV{CC}")
else()
    message("Environment variables cc and/or CC are not set, using default settings")
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

file(GLOB SOURCEFILES "src/*.cpp")
include_directories(BEFORE PRIVATE "include" "../include" "nlopt/include" "nlopt/build" "/usr/include/c++/12")
link_directories(BEFORE  PRIVATE "nlopt/build")
link_libraries(nlopt m)
add_executable(nwq_vqe "main.cpp" ${SOURCEFILES})
add_library(vqe SHARED ${SOURCEFILES})



add_subdirectory("examples")

