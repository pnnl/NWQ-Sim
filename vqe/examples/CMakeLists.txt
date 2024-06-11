cmake_minimum_required(VERSION 3.20)


include_directories(BEFORE PRIVATE "include" "../include" ${NLOPT_DIR})
link_directories(BEFORE  PRIVATE ${NLOPT_DIR})
link_libraries(vqe)

add_executable(basic_example_cuda "basic_example_cuda.cu" ${SOURCE})
add_executable(basic_example_cpu "basic_example_cpu.cpp" ${SOURCE})
add_executable(config_example "config_example.cpp" ${SOURCE})
