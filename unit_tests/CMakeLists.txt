#set(GTEST_INCLUDE_DIR /usr/local/include/google)
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIR})

add_subdirectory(gtest-1.7.0)
include_directories(gtest-1.7.0/include)
set(GTEST_BOTH_LIBRARIES gtest gtest_main)

add_subdirectory(simulation)   
add_subdirectory(math)   
add_subdirectory(calib) 
add_subdirectory(io)  
