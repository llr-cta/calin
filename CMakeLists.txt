cmake_minimum_required(VERSION 3.0)  
project(calin)
enable_testing()

# ======================================
# Do some annoying incomprehensible CMake crap
# ======================================

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/CMakeModules/")

include_directories(include ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)

## include FFTW libraries
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${FFTW_LIBRARIES} ) 

## include GSL libraries
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${GSL_LIBRARIES} ) 

add_subdirectory(proto)
add_subdirectory(src)
add_subdirectory(unit_tests)
