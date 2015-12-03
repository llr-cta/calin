cmake_minimum_required(VERSION 3.0)  
project(calin)
enable_testing()

# ======================================
# Do some annoying incomprehensible CMake crap
# ======================================

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/CMakeModules/")

include_directories(include ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

### include Protobuf libraries
find_package(Protobuf REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIRS} ${CMAKE_CURRENT_BUILD_DIR}) 
set(EXTRA_LIBS ${EXTRA_LIBS} ${PROTOBUF_LIBRARY})

## include FFTW libraries
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${FFTW_LIBRARIES} ) 

## include GSL libraries
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${GSL_LIBRARIES} ) 

## include Sqlite3 libraries
find_package(SQLite3 REQUIRED)
include_directories(${SQLITE3_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${SQLITE3_LIBRARY} ) 

#Path to a python
#set(PYTHON_EXECUTABLE /opt/local/bin/python3.4)
#set(PYTHON_INCLUDE_DIR /opt/local/Library/Frameworks/Python.framework/Versions/3.4/include/python3.4m)
#set(PYTHON_LIBRARY /opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/config-3.4m/libpython3.4.dylib)

find_package(PythonInterp) 
find_package(PythonLibs) 

## Find Geant4
if(NO_GEANT4)
	message(STATUS "Skipping GEANT4 as requested")
	SET(GEANT4_FOUND FALSE)
else()
       find_package(Geant4) 
       if(GEANT4_FOUND)
		message(STATUS "GEANT4_INCLUDE_DIR: " ${GEANT4_INCLUDE_DIR}) 
		message(STATUS "GEANT4_LIBRARY_DIR: " ${GEANT4_LIBRARY_DIR}) 
		include_directories(${GEANT4_INCLUDE_DIR})
		link_directories(${GEANT4_LIBRARY_DIR}) 
	endif()
endif()

add_subdirectory(tools)
add_subdirectory(proto)
add_subdirectory(src)  
add_subdirectory(swig)  
add_subdirectory(unit_tests)
