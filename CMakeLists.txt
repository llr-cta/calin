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

#Path to a program.
set(PYTHON_EXECUTABLE /opt/local/bin/python3.4)

#Path to a file.
set(PYTHON_INCLUDE_DIR /opt/local/Library/Frameworks/Python.framework/Versions/3.4/include/python3.4m)

#Path to a library.
set(PYTHON_LIBRARY /opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/config-3.4m/libpython3.4.dylib)

find_package(PythonLibs)

add_subdirectory(proto)
add_subdirectory(src)  
add_subdirectory(swig)  
add_subdirectory(unit_tests)
