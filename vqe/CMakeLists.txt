cmake_minimum_required(VERSION 3.29.3)
project(nwq_vqe)


set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 17)
if (NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

file(GLOB SOURCEFILES "src/*.cpp")
add_executable(nwq_vqe ${SOURCEFILES})

target_include_directories(nwq_vqe BEFORE PRIVATE "include" "../include" ${NLOPT_DIR})
target_link_directories(nwq_vqe BEFORE  PRIVATE ${NLOPT_DIR})
target_link_libraries(nwq_vqe nlopt m)

