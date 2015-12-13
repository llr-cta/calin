PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 calin.proto common_types.proto unittest.proto)

message(STATUS "PROTOBUF_IMPORT_DIRS: " ${PROTOBUF_IMPORT_DIRS}) 
message(STATUS "PROTOBUF_INCLUDE_DIRS: " ${PROTOBUF_INCLUDE_DIRS}) 

set(CALIN_TARGET_LIBRARY calin_proto)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${ProtoSources}) 
include_directories(${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CALIN_TARGET_LIBRARY} ${PROTOBUF_LIBRARY}) 
install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)

add_subdirectory(math)
add_subdirectory(simulation)
add_subdirectory(io)

