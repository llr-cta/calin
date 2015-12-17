include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CALIN_MATH_SWIG_FILES function.i pdf_1d.i histogram.i optimizer.i
			  hessian.i special.i hex_array.i rng.i vs_physics.i)

set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES CPLUSPLUS ON)
#set_source_files_properties(${CALIN_MATH_SWIG_FILES} PROPERTIES SWIG_FLAGS "-importall")

foreach (interface_file ${CALIN_MATH_SWIG_FILES})
	string(REPLACE ".i" "" output_file ${interface_file})
#	message(STATUS ${interface_file} ${output_file}) 
	swig_add_module(${output_file} python ${interface_file})
	swig_link_libraries(${output_file} calin_math ${PYTHON_LIBRARIES})
	install(TARGETS ${SWIG_MODULE_${output_file}_REAL_NAME}  DESTINATION python/calin/math)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_file}.py  DESTINATION python/calin/math)
endforeach (interface_file)

add_dependencies(_histogram _histogram_pb)

install(FILES __init__.py DESTINATION python/calin/math)
