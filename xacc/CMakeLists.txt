find_package(XACC REQUIRED)

add_executable(xacc-bell bell.cpp)

target_link_libraries(xacc-bell ${XACC_LIBRARIES})