#add_subdirectory(math)
#add_subdirectory(common)

set(PROTOBUF_IMPORT_DIRS /usr/local/include ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Protobuf REQUIRED)

PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders
					 cta.proto histogram.proto)

set(PROTOBUF_INCLUDE_DIRS /usr/local/include) 

#add_library(calin_proto STATIC)
add_library(calin_proto STATIC ${ProtoSources}) 

include_directories(${PROTOBUF_INCLUDE_DIRS} ${CMAKE_CURRENT_BUILD_DIR}) 

target_link_libraries(calin_proto ${PROTOBUF_LIBRARY}) 

