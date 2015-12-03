add_executable(protoc-gen-swig protoc-gen-swig.cpp swig_generator.cpp)

include_directories(${CMAKE_CURRENT_BUILD_DIR}) 

target_link_libraries(protoc-gen-swig ${PROTOBUF_LIBRARY}
				${PROTOBUF_PROTOC_LIBRARIES}) 

install(TARGETS protoc-gen-swig DESTINATION bin)
