find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

find_package(Numpy)
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

set(CMAKE_SWIG_FLAGS "")
add_compile_options(-DSWIG_TYPE_TABLE=calin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}) 

#swig_add_moduleSWIG_ADD_MODULE(example python example.i example.cxx)
#SWIG_LINK_LIBRARIES(example ${PYTHON_LIBRARIES})

add_subdirectory(math)
