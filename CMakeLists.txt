
#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project(EnGPar)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

#Setup doxygen documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
    )
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif()

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

option(ENABLE_PUMI "Build with PUMI interface" OFF)

if (ENABLE_PUMI)
  find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
  find_package(SCOREC 2.1.0 REQUIRED COMPONENTS pcu CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
endif()

option(ENABLE_ZOLTAN "Enable Zoltan iterface [ON|OFF]" OFF)
message(STATUS "ENABLE_ZOLTAN: " ${ENABLE_ZOLTAN})
if (ENABLE_ZOLTAN)
  find_package(Zoltan MODULE REQUIRED)
endif()

option(ENABLE_PARMETIS "Enable ParMETIS iterface [ON|OFF]" OFF)
message(STATUS "ENABLE_PARMETIS: " ${ENABLE_PARMETIS})
if (ENABLE_PARMETIS)
  find_package(Parmetis MODULE REQUIRED)
  add_definitions(-DHAS_PARMETIS)
endif()

option(ENGPAR_FORTRAN_INTERFACE "enable fortran interface" OFF)
message(STATUS "ENGPAR_FORTRAN_INTERFACE: ${ENGPAR_FORTRAN_INTERFACE}")
if(ENGPAR_FORTRAN_INTERFACE)
  enable_language(Fortran)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
endif()

if (ENABLE_KOKKOS)
  find_package(Kokkos MODULE REQUIRED)
  add_definitions(-DKOKKOS_ENABLED)
endif()

include_directories(
  ${SCOREC_PREFIX}/include
  ${KOKKOS_INCLUDE_DIRS}
  )

macro(engpar_export_lib target headers)
  install(TARGETS ${target}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
  install(FILES ${headers} DESTINATION include)
endmacro(engpar_export_lib)

add_subdirectory(support)

add_subdirectory(agi)

add_subdirectory(ssg)

add_subdirectory(zoltan)

add_subdirectory(partition)

add_subdirectory(interfaces)

if(ENGPAR_FORTRAN_INTERFACE)
  add_subdirectory(ftn)
endif()

add_subdirectory(test)
