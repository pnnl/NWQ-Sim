set(PYBIND11_FIND_STRATEGY "CONFIG")

if(NOT TARGET pybind11::module)
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        if(EXISTS "${PROJECT_SOURCE_DIR}/vqe/pybind11/CMakeLists.txt")
            add_subdirectory(${PROJECT_SOURCE_DIR}/vqe/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/pybind11)
        else()
            message(FATAL_ERROR "pybind11 not found. Please install pybind11 or initialize the vqe/pybind11 submodule.")
        endif()
    endif()
endif()

pybind11_add_module(nwqsim MODULE src/nwqsim_py.cpp)

target_compile_features(nwqsim PRIVATE cxx_std_17)

target_include_directories(nwqsim PRIVATE ${PROJECT_SOURCE_DIR}/include)

if(OpenMP_FOUND)
    target_link_libraries(nwqsim PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(nwqsim PRIVATE OMP_ENABLED)
endif()

if(MPI_FOUND)
    target_link_libraries(nwqsim PRIVATE MPI::MPI_CXX)
    target_compile_definitions(nwqsim PRIVATE MPI_ENABLED)
endif()

if(CUDAToolkit_FOUND)
    target_compile_definitions(nwqsim PRIVATE CUDA_ENABLED)
endif()

if(HIP_FOUND)
    target_compile_definitions(nwqsim PRIVATE HIP_ENABLED)
endif()

set_target_properties(nwqsim PROPERTIES OUTPUT_NAME "nwqsim")
