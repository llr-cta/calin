#add_compile_options(-DSWIG_TYPE_TABLE=calin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}) 

set(CALIN_SWIG_FILES package_wide_definitions.i google_protobuf.i)

set_source_files_properties(${CALIN_SWIG_FILES} PROPERTIES CPLUSPLUS ON)

foreach (interface_file ${CALIN_SWIG_FILES})
        string(REPLACE ".i" "" output_file ${interface_file})
#       message(STATUS ${interface_file} ${output_file})
        swig_add_module(${output_file} python ${interface_file})
        swig_link_libraries(${output_file} ${PYTHON_LIBRARIES} ${PROTOBUF_LIBRARIES})
        install(TARGETS ${SWIG_MODULE_${output_file}_REAL_NAME}  DESTINATION python/calin)
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${output_file}.py  DESTINATION python/calin)
endforeach (interface_file)

install(FILES __init__.py DESTINATION python/calin)

add_subdirectory(math) 
add_subdirectory(io) 
add_subdirectory(calib) 
add_subdirectory(simulation) 

install(FILES __init__.py DESTINATION python/calin)
