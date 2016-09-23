project (TestSuite CXX)
include_directories (${CMAKE_SOURCE_DIR}/include/)
add_executable (tully_1.exe ./tully_1/main.cpp)
target_link_libraries (tully_1.exe scattering)

add_custom_target (test_suite)
add_dependencies (test_suite tully_1.exe)
