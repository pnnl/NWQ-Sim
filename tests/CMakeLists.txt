enable_testing()

# Assume that we have a test file test1.cpp

add_executable(circuit_test circuit_test.cpp)
add_executable(cpu_sin_test cpu_sin_test.cpp)

# Include directories for test1 target

target_include_directories(circuit_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(cpu_sin_test PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(OpenMP_FOUND)
    target_link_libraries(circuit_test PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(cpu_sin_test PUBLIC OpenMP::OpenMP_CXX)
endif()

# Add the test
add_test(NAME Test1 COMMAND circuit_test)
add_test(NAME Test2 COMMAND cpu_sin_test)