
#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project(EnGPar)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

#Settings options for testing
#enable_testing()
#include(CTest)
#TODO: Locate scorec directory. Module? FindSCOREC.cmake? or whatever
find_package(SCOREC MODULE REQUIRED)

include_directories(${SCOREC_INCLUDE_DIRS})

add_subdirectory(agi)

add_subdirectory(interfaces)

add_subdirectory(test)
