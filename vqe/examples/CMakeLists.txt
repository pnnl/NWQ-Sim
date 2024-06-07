cmake_minimum_required(VERSION 3.20)

file(GLOB SOURCE "../src/*.cpp")

include_directories(BEFORE PRIVATE "include" "../include" ${NLOPT_DIR})
link_directories(BEFORE  PRIVATE ${NLOPT_DIR})
link_libraries(vqe)

add_executable(basic_example "basic_example.cpp" ${SOURCE})
