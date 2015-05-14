include_directories(${PYTHON_INCLUDE_PATH})
add_library(calin_io SHARED log.cpp)
target_link_libraries(calin_io ${PYTHON_LIBRARIES})
install(TARGETS calin_io DESTINATION lib)
