cmake_minimum_required (VERSION 3.5)

project (LibScattering CXX)

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

include (BasicSetup)

set (${PROJNAME}_EIGEN_INCLUDE /usr/local/include/ CACHE PATH
  "path to the Eigen headers")

set (SRC_FILES ./src/scattering_1d.cpp ./src/scattering.cpp)

add_library (scattering SHARED ${SRC_FILES})
target_include_directories (scattering PRIVATE ./include)
set_target_properties (scattering PROPERTIES 
  CXX_STANDARD ${${PROJNAME}_CXXSTANDARD}
  CXX_STANDARD_REQUIRED TRUE
  CXX_VISIBILITY_PRESET hidden
  VERSION 1.0.0
  SOVERSION 1.0
  )
target_compile_options (scattering PUBLIC ${${PROJNAME}_WARNING_FLAG})
target_include_directories (scattering PUBLIC ${${PROJNAME}_EIGEN_INCLUDE})

add_library (scattering_static STATIC ${SRC_FILES})
target_include_directories (scattering_static PRIVATE ./include)
set_target_properties (scattering_static PROPERTIES 
  CXX_STANDARD ${${PROJNAME}_CXXSTANDARD}
  CXX_STANDARD_REQUIRED TRUE
  CXX_VISIBILITY_PRESET hidden
  )
target_compile_options (scattering_static PUBLIC ${${PROJNAME}_WARNING_FLAG})
target_include_directories (scattering_static
  PUBLIC ${${PROJNAME}_EIGEN_INCLUDE})

install (TARGETS scattering scattering_static LIBRARY DESTINATION lib/
  ARCHIVE DESTINATION lib/)
install (FILES ./include/scattering_1d.hpp DESTINATION include/)


add_subdirectory (examples EXCLUDE_FROM_ALL)
