cmake_minimum_required(VERSION 3.20)

project(nwqsim_refactored LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(nwqsim_refactored STATIC
    simulator/fusion.cpp
    simulator/sv_sim.cpp
)

target_include_directories(nwqsim_refactored
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

add_executable(refactored_statevector_test
    test_cpu_state.cpp
)

target_link_libraries(refactored_statevector_test
    PRIVATE
        nwqsim_refactored
)

target_include_directories(refactored_statevector_test
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)
