
#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project(EnGPar)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")


#Settings options for testing
enable_testing()
include(CTest)

option(IS_TESTING "Build for CTest" OFF)
message(STATUS "IS_TESTING: ${IS_TESTING}")

get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH)
find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
set(MPIRUN_PROCFLAG "-np" CACHE string
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")
set(VALGRIND "" CACHE string
  "the valgrind executable")
set(VALGRIND_ARGS "" CACHE string
  "the command line arguments to VALGRIND")

set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")
find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)

include_directories(${SCOREC_PREFIX}/include)

add_subdirectory(agi)

add_subdirectory(zoltan)

add_subdirectory(interfaces)

add_subdirectory(test)
