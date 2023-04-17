# - Try to find FFTW3.

# set the library
set(_libraries)
list(APPEND _libraries fftw3)

# Keep a list of variable names 
set(_check_list)

# Search for all requested libraries.
foreach(_lib ${_libraries})
  string(TOUPPER ${_lib} _LIB)
  find_library(${_LIB}_LIBRARY ${_lib}
    HINTS ${FFTW3_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(${_LIB}_LIBRARY)
  list(APPEND FFTW3_LIBRARIES ${${_LIB}_LIBRARY})
  list(APPEND _check_list ${_LIB}_LIBRARY)
endforeach(_lib ${_libraries})

# Search for the header file.
find_path(FFTW3_INCLUDE_DIR fftw3.h 
  HINTS ${FFTW3_ROOT_DIR} PATH_SUFFIXES include)
mark_as_advanced(FFTW3_INCLUDE_DIR)
list(APPEND _check_list FFTW3_INCLUDE_DIR)

# Handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 DEFAULT_MSG ${_check_list})
