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
  VERSION 1.0.0
  SOVERSION 1.0
  )
target_compile_options (scattering PUBLIC ${${PROJNAME}_WARNING_FLAG})
target_link_libraries (scattering ${${PROJNAME}_MKL_LIBS})

install (TARGETS scattering LIBRARY DESTINATION lib/)
install (FILES ./include/scattering_1d.hpp DESTINATION include/)

add_subdirectory (examples EXCLUDE_FROM_ALL)

# check MKL header
include (CheckCXXSourceCompiles)
set (CMAKE_REQUIRED_INCLUDES ${${PROJNAME}_MKL_INCLUDE})
unset (IS_MKL_H_AVAILABLE CACHE)
CHECK_CXX_SOURCE_COMPILES ("#include <mkl.h>\nint main() { return 0; }"
  IS_MKL_H_AVAILABLE)
if (NOT IS_MKL_H_AVAILABLE)
  message (WARNING "mkl.h is not found. "
    "Please change ${PROJNAME}_MKL_INCLUDE to point to the header directory.")
endif ()

# check PARDISO
if (IS_MKL_H_AVAILABLE)
  include (CheckCXXSymbolExists)
  set (CMAKE_REQUIRED_LIBRARIES ${${PROJNAME}_MKL_LIBS})
  set (CMAKE_REQUIRED_FLAGS ${${PROJNAME}_MKL_LINKFLAG})
  unset (IS_PARDISO_AVAILABLE CACHE)
  CHECK_CXX_SYMBOL_EXISTS ("pardiso" "mkl.h" IS_PARDISO_AVAILABLE)
  if (NOT IS_PARDISO_AVAILABLE)
    message (WARNING "PARDISO does not seem to be available. "
      "Your MKL might be too old.")
  endif ()
endif ()
