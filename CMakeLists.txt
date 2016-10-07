cmake_minimum_required (VERSION 3.5)

project (LibScattering CXX)

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

include (BasicSetup)

set (${PROJNAME}_MKL_LINKFLAG "-L$ENV{MKLROOT}/lib/intel64"
  CACHE STRING "link flag for Intel MKL")
set (${PROJNAME}_MKL_LIBS "mkl_intel_lp64;mkl_core;mkl_sequential;pthread;m;dl"
  CACHE STRING "Intel MKL libraries")
set (${PROJNAME}_MKL_INCLUDE "$ENV{MKLROOT}/include"
  CACHE PATH "Intel MKL header location")

add_library (scattering SHARED ./src/scattering_1d.cpp
                                  ./src/scattering.cpp)
target_include_directories (scattering PRIVATE ./include
  PRIVATE ${${PROJNAME}_MKL_INCLUDE})
set_target_properties (scattering PROPERTIES
  LINK_FLAGS ${${PROJNAME}_MKL_LINKFLAG}
  CXX_STANDARD ${${PROJNAME}_CXXSTANDARD}
  CXX_STANDARD_REQUIRED TRUE
  CXX_VISIBILITY_PRESET hidden
  )
target_compile_options (scattering PUBLIC ${${PROJNAME}_WARNING_FLAG})
target_link_libraries (scattering ${${PROJNAME}_MKL_LIBS})


add_subdirectory (examples EXCLUDE_FROM_ALL)
