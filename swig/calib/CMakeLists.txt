include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CALIN_SWIG_FILES spe_fit.i)
set(CALIN_SWIG_TARGET calib)

set_source_files_properties(${CALIN_SWIG_FILES} PROPERTIES CPLUSPLUS ON)
#set_source_files_properties(${CALIN_SWIG_FILES} PROPERTIES SWIG_FLAGS "-importall")

foreach (interface_file ${CALIN_SWIG_FILES})
	string(REPLACE ".i" "" output_file ${interface_file})
#	message(STATUS ${interface_file} ${output_file})
	swig_add_module(${output_file} python ${interface_file})
	swig_link_libraries(${output_file} calin_${CALIN_SWIG_TARGET} ${PYTHON_LIBRARIES})
	install(TARGETS ${SWIG_MODULE_${output_file}_REAL_NAME}  DESTINATION python/calin/${CALIN_SWIG_TARGET})
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_file}.py  DESTINATION python/calin/${CALIN_SWIG_TARGET})
endforeach (interface_file)

install(FILES __init__.py DESTINATION python/calin/${CALIN_SWIG_TARGET})
