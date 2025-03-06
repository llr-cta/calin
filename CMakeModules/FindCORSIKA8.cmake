# - Find CORSIKA8
# Find the native CORSIKA8 includes and library
# This module defines
#  CORSIKA8_INCLUDE_DIR, where to find fftw3.h, etc.
#  CORSIKA8_LIBRARIES, the libraries needed to use CORSIKA8.
#  CORSIKA8_FOUND, If false, do not try to use CORSIKA8.
# also defined, but not for general use are
#  CORSIKA8_LIBRARY, where to find the CORSIKA8 library.


find_path ( CORSIKA8_INCLUDE_DIR corsika.hpp PATH_SUFFIXES corsika HINTS ${CORSIKA8_DIR} )
cmake_path( GET CORSIKA8_INCLUDE_DIR PARENT_PATH CORSIKA8_INCLUDE_DIR)

find_library ( CORSIKA8_DATA_LIBRARY PATH_SUFFIXES corsika NAMES CorsikaData HINTS ${CORSIKA8_DIR} )
find_library ( PROPOSAL_LIBRARY NAMES PROPOSAL HINTS ${CORSIKA8_DIR} )
find_library ( EPOS_LIBRARY PATH_SUFFIXES corsika NAMES EPOS_static HINTS ${CORSIKA8_DIR} )
find_library ( QGSJET_LIBRARY PATH_SUFFIXES corsika NAMES QGSJetII_static HINTS ${CORSIKA8_DIR} )
find_library ( SIBYLL_LIBRARY PATH_SUFFIXES corsika NAMES Sibyll HINTS ${CORSIKA8_DIR} )
find_library ( PYTHIA8_LIBRARY PATH_SUFFIXES corsika NAMES pythia8 HINTS ${CORSIKA8_DIR} )
find_library ( CUBIC_INTERPOLATION_LIBRARY NAMES CubicInterpolation HINTS ${CORSIKA8_DIR} )
find_library ( SPDLOG_LIBRARY NAMES spdlog HINTS ${CORSIKA8_DIR} )
find_library ( FMT_LIBRARY NAMES fmt HINTS ${CORSIKA8_DIR} )
find_library ( YAML_LIBRARY NAMES yaml-cpp HINTS ${CORSIKA8_DIR} )
find_library ( BZ2_LIBRARY NAMES bz2 HINTS ${CORSIKA8_DIR} )
find_library ( GFORTRAN_LIBRARY NAMES gfortran HINTS ${CORSIKA8_DIR} /opt/local/lib/gcc13 )

set ( CORSIKA8_LIBRARIES 
        ${CORSIKA8_DATA_LIBRARY} 
        ${PROPOSAL_LIBRARY} 
        ${EPOS_LIBRARY}
        ${QGSJET_LIBRARY}
        ${SIBYLL_LIBRARY}
        ${PYTHIA8_LIBRARY}
        ${CUBIC_INTERPOLATION_LIBRARY} 
        ${SPDLOG_LIBRARY}
        ${FMT_LIBRARY}
        ${YAML_LIBRARY}
        ${BZ2_LIBRARY}
        ${GFORTRAN_LIBRARY} )
set ( CORSIKA8_INCLUDE_DIRS ${CORSIKA8_INCLUDE_DIR} )

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set CORSIKA8_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( CORSIKA8 DEFAULT_MSG 
        CORSIKA8_DATA_LIBRARY 
        PROPOSAL_LIBRARY 
        EPOS_LIBRARY 
        QGSJET_LIBRARY
        SIBYLL_LIBRARY
        PYTHIA8_LIBRARY
        CUBIC_INTERPOLATION_LIBRARY 
        SPDLOG_LIBRARY 
        FMT_LIBRARY 
        YAML_LIBRARY 
        BZ2_LIBRARY
        GFORTRAN_LIBRARY
        CORSIKA8_INCLUDE_DIR )
