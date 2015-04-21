include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CALIN_MATH_SWIG_FILES function.i pdf_1d.i)

set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES CPLUSPLUS ON)
#set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES SWIG_FLAGS "-includeall")

foreach (interface_file ${CALIN_MATH_SWIG_FILES})
	string(REPLACE ".i" "" output_file ${interface_file})
#	message(STATUS ${interface_file} ${output_file}) 
	swig_add_module(${output_file} python ${interface_file})
	swig_link_libraries(${output_file} calin_math ${PYTHON_LIBRARIES})
endforeach (interface_file)

