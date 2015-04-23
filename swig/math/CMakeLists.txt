include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CALIN_MATH_SWIG_FILES function.i pdf_1d.i histogram.i optimizer.i)

set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES CPLUSPLUS ON)
#set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES SWIG_FLAGS "-py3")

foreach (interface_file ${CALIN_MATH_SWIG_FILES})
	string(REPLACE ".i" "" output_file ${interface_file})
#	message(STATUS ${interface_file} ${output_file}) 
	swig_add_module(${output_file} python ${interface_file})
	swig_link_libraries(${output_file} calin_math ${PYTHON_LIBRARIES})
	install(TARGETS ${SWIG_MODULE_${output_file}_REAL_NAME}  DESTINATION python/calin/math)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_file}.py  DESTINATION python/calin/math)
endforeach (interface_file)

install(FILES __init__.py DESTINATION python/calin/math)
