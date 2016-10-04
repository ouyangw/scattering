project (Examples CXX)
include_directories (${CMAKE_SOURCE_DIR}/include/)

add_executable (tully_1.exe ./tully_1/main.cpp)
target_link_libraries (tully_1.exe scattering)
set_target_properties (tully_1.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/)

add_executable (tully_1_range.exe ./tully_1_range/main.cpp)
target_link_libraries (tully_1_range.exe scattering)
set_target_properties (tully_1_range.exe PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY bin/)

add_executable (tully_1_print.exe ./tully_1_print/main.cpp)
target_link_libraries (tully_1_print.exe scattering)
set_target_properties (tully_1_print.exe PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY bin/)

add_executable (tully_2.exe ./tully_2/main.cpp)
target_link_libraries (tully_2.exe scattering)
set_target_properties (tully_2.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/)

add_executable (tully_2_range.exe ./tully_2_range/main.cpp)
target_link_libraries (tully_2_range.exe scattering)
set_target_properties (tully_2_range.exe PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY bin/)

add_executable (tully_3.exe ./tully_3/main.cpp)
target_link_libraries (tully_3.exe scattering)
set_target_properties (tully_3.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/)

add_executable (tully_3_range.exe ./tully_3_range/main.cpp)
target_link_libraries (tully_3_range.exe scattering)
set_target_properties (tully_3_range.exe PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY bin/)

add_custom_target (examples)
add_dependencies (examples tully_1.exe
                           tully_1_range.exe
                           tully_1_print.exe
                           tully_2.exe
                           tully_2_range.exe
                           tully_3.exe
                           tully_3_range.exe
                           )
