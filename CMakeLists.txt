cmake_minimum_required (VERSION 3.5)

project (LibScattering CXX)

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

include (BasicSetup)

set (${PROJNAME}_EIGEN_INCLUDE /usr/local/include/ CACHE PATH
  "path to the Eigen headers")

add_library (scattering SHARED ./src/scattering_1d.cpp
                                  ./src/scattering.cpp)
target_include_directories (scattering PRIVATE ./include)
set_target_properties (scattering PROPERTIES 
  CXX_STANDARD ${${PROJNAME}_CXXSTANDARD}
  CXX_STANDARD_REQUIRED TRUE
  CXX_VISIBILITY_PRESET hidden
  )
target_compile_options (scattering PUBLIC ${${PROJNAME}_WARNING_FLAG})
target_include_directories (scattering PUBLIC ${${PROJNAME}_EIGEN_INCLUDE})


add_subdirectory (examples EXCLUDE_FROM_ALL)
