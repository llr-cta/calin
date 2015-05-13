add_library(calin_io SHARED log.cpp)
target_link_libraries(calin_io)

install(TARGETS calin_io DESTINATION lib)
