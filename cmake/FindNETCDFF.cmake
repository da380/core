# - Try to find NETCDFF.

# Variables used by this module:
#  NETCDFF_ROOT_DIR             - NETCDFF root directory
# Variables defined by this module:
#  NETCDFF_FOUND                - system has NETCDFF
#  NETCDFF_INCLUDE_DIR          - the NETCDFF include directory (cached)
#  NETCDFF_INCLUDE_DIRS         - the NETCDFF include directories
#                               (identical to NETCDFF_INCLUDE_DIR)
#  NETCDFF_LIBRARY              - the NETCDFF library (cached)
#  NETCDFF_LIBRARIES            - list of all NETCDFF libraries found




# Loop over each component.
set(_libraries)
list(APPEND _libraries netcdff)

# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
set(_check_list)

# Search for all requested libraries.
foreach(_lib ${_libraries})
  string(TOUPPER ${_lib} _LIB)
  find_library(${_LIB}_LIBRARY ${_lib}
    HINTS ${NETCDFF_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(${_LIB}_LIBRARY)
  list(APPEND NETCDFF_LIBRARIES ${${_LIB}_LIBRARY})
  list(APPEND _check_list ${_LIB}_LIBRARY)
endforeach(_lib ${_libraries})

# Search for the header file.
find_path(NETCDFF_INCLUDE_DIR netcdf.mod 
  HINTS ${NETCDFF_ROOT_DIR} PATH_SUFFIXES include)
mark_as_advanced(NETCDFF_INCLUDE_DIR)
list(APPEND _check_list NETCDFF_INCLUDE_DIR)
set(NETCDFF_MODULE NETCDFF_INCLUDE_DIR/netcdf.mod)


# Handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NETCDFF DEFAULT_MSG ${_check_list})
