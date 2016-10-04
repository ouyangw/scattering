cmake_minimum_required (VERSION 3.5)

project (LibScattering CXX)

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

include (BasicSetup)

set (${PROJNAME}_MKL_LDFLAG "-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl" CACHE STRING "link flag for Intel MKL")

add_library (scattering SHARED ./src/scattering_1d.cpp
                                  ./src/scattering.cpp)
target_include_directories (scattering PRIVATE ./include)
set_target_properties (scattering PROPERTIES 
  LINK_FLAGS ${${PROJNAME}_MKL_LDFLAG}
  CXX_STANDARD ${${PROJNAME}_CXXSTANDARD}
  CXX_STANDARD_REQUIRED TRUE
  CXX_VISIBILITY_PRESET hidden
  )
target_compile_options (scattering PUBLIC ${${PROJNAME}_WARNING_FLAG})


add_subdirectory (examples EXCLUDE_FROM_ALL)
