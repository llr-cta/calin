find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

#Path to a program.
set(PYTHON_EXECUTABLE /opt/local/bin/python3.4)

#Path to a file.
set(PYTHON_INCLUDE_DIR /opt/local/Library/Frameworks/Python.framework/Versions/3.4/include/python3.4m)

#Path to a library.
set(PYTHON_LIBRARY /opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/config-3.4m/libpython3.4.dylib)


find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

find_package(Numpy)
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

set(CMAKE_SWIG_FLAGS "")
#add_compile_options(-DSWIG_TYPE_TABLE=calin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}) 

#swig_add_moduleSWIG_ADD_MODULE(example python example.i example.cxx)
#SWIG_LINK_LIBRARIES(example ${PYTHON_LIBRARIES})

add_subdirectory(math) 
add_subdirectory(calib) 

install(FILES __init__.py DESTINATION python/calin)
