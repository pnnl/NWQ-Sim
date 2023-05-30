enable_testing()

# Assume that we have a test file test1.cpp

add_executable(circuit_test circuit_test.cpp)
add_executable(cpu_sin cpu_sin.cpp)

# Include directories for test1 target

target_include_directories(circuit_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(cpu_sin PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(OpenMP_FOUND)
    target_link_libraries(circuit_test PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(cpu_sin PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MPI_FOUND)
    add_executable(cpu_mpi cpu_mpi.cpp)
    target_include_directories(cpu_mpi PRIVATE ${MPI_INCLUDE_PATH})
    target_link_libraries(cpu_mpi PUBLIC MPI::MPI_CXX)
endif()

# Add the test
# add_test(NAME Test1 COMMAND circuit_test)
# add_test(NAME Test2 COMMAND cpu_sin)