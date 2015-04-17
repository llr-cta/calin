include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CALIN_MATH_SWIG_FILES function.i pdf_1d.i)

set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES CPLUSPLUS ON)
#set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES SWIG_FLAGS "-includeall")

swig_add_module(function python function.i)
swig_link_libraries(function calin_math ${PYTHON_LIBRARIES})

