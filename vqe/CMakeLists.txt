cmake_minimum_required(VERSION 3.29.3)
project(nwq_vqe)


set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 17)
if (NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

file(GLOB SOURCEFILES "src/*.cpp")
include_directories(BEFORE PRIVATE "include" "../include" ${NLOPT_DIR})
link_directories(BEFORE  PRIVATE ${NLOPT_DIR})
link_libraries(nlopt m)
add_executable(nwq_vqe "main.cpp" ${SOURCEFILES})
add_library(vqe SHARED ${SOURCEFILES})



add_subdirectory("examples")

