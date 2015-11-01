#add_subdirectory(math)
#add_subdirectory(common)

find_package(Protobuf REQUIRED)

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 calin.proto calin_common_types.proto
					 io_sql_transceiver.proto
					 math_histogram.proto unittest.proto)

message(STATUS "PROTOBUF_IMPORT_DIRS: " ${PROTOBUF_IMPORT_DIRS}) 
message(STATUS "PROTOBUF_INCLUDE_DIRS: " ${PROTOBUF_INCLUDE_DIRS}) 

add_library(calin_proto SHARED ${ProtoSources}) 

include_directories(${PROTOBUF_INCLUDE_DIRS} ${CMAKE_CURRENT_BUILD_DIR}) 

target_link_libraries(calin_proto ${PROTOBUF_LIBRARY}) 

install(TARGETS calin_proto DESTINATION lib)

