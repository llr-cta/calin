# - Find CORSIKA8
# Find the native CORSIKA8 includes and library
# This module defines
#  CORSIKA8_INCLUDE_DIR, where to find fftw3.h, etc.
#  CORSIKA8_LIBRARIES, the libraries needed to use CORSIKA8.
#  CORSIKA8_FOUND, If false, do not try to use CORSIKA8.
# also defined, but not for general use are
#  CORSIKA8_LIBRARY, where to find the CORSIKA8 library.

find_path ( CORSIKA8_INCLUDE_DIR corsika.hpp PATH_SUFFIXES corsika )
find_library ( CORSIKA8_DATA_LIBRARY PATH_SUFFIXES corsika NAMES CorsikaData  )

set ( CORSIKA8_LIBRARIES ${CORSIKA8_DATA_LIBRARY}  )
set ( CORSIKA8_INCLUDE_DIRS ${CORSIKA8_INCLUDE_DIR} )

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set CORSIKA8_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( CORSIKA8 DEFAULT_MSG CORSIKA8_DATA_LIBRARY CORSIKA8_INCLUDE_DIR )
