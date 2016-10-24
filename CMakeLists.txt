
#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project(EnGPar)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

#Settings options for testing
#enable_testing()
#include(CTest)
set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")
find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)

include_directories(${SCOREC_PREFIX}/include)

add_subdirectory(agi)

add_subdirectory(interfaces)

add_subdirectory(test)
