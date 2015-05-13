set(GTEST_INCLUDE_DIR /usr/local/include/google)
find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIR})

add_subdirectory(math)  
add_subdirectory(calib) 
add_subdirectory(io)  
