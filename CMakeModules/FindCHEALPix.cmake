# - Find CHEALPIX
# Find the native CHEALPix includes and library
# This module defines
#  CHEALPIX__INCLUDE_DIR, where to find chealpix.h, etc.
#  CHEALPIX_LIBRARIES, the libraries needed to use CHEALPix.
#  CHEALPIX_FOUND, If false, do not try to use CHEALPix.
# also defined, but not for general use are
#  CHEALPIX_LIBRARY, where to find the CHEALPix library.

find_path ( CHEALPIX_INCLUDE_DIR chealpix.h )
find_library ( CHEALPIX_LIBRARY NAMES chealpix )

set ( CHEALPIX_LIBRARIES ${CHEALPIX_LIBRARY} )
set ( CHEALPIX_INCLUDE_DIRS ${CHEALPIX_INCLUDE_DIR} )

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( CHEALPIX DEFAULT_MSG CHEALPIX_LIBRARY CHEALPIX_INCLUDE_DIR )
