#Setting Version Number, Project Name
cmake_minimum_required (VERSION 3.6)
project(EnGPar VERSION 1.1.0 LANGUAGES CXX C)

include(cmake/bob.cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

bob_begin_package()

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

execute_process(COMMAND ${CMAKE_CXX_COMPILER} -show OUTPUT_VARIABLE
  VERBOSE_CXX_COMPILE_CMD)
message(STATUS "CXX compilation command: ${VERBOSE_CXX_COMPILE_CMD}")
get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH)
find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
set(MPIRUN_PROCFLAG "-np" CACHE STRING
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")
set(VALGRIND "" CACHE STRING
  "the valgrind executable")
set(VALGRIND_ARGS "" CACHE STRING
  "the command line arguments to VALGRIND")

if (ENABLE_KOKKOS)
  find_package(KokkosCore REQUIRED)
  find_package(KokkosKernels REQUIRED)
  add_definitions(-DENGPAR_KOKKOS_ENABLED)
  set(ENGPAR_KOKKOS_ENABLED true)
  include_directories(${KokkosCore_INCLUDE_DIRS})

  message(STATUS "Kokkos enabled - Overriding compiler selection for consistency with Kokkos install")
  set(CMAKE_CXX_COMPILER ${KokkosCore_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${KokkosCore_C_COMPILER})
  set(CMAKE_FORTRAN_COMPILER ${KokkosCore_FORTRAN_COMPILER})
  message(STATUS "Kokkos CMAKE_CXX_COMPILER ${KokkosCore_CXX_COMPILER}")
  message(STATUS "Kokkos CMAKE_C_COMPILER ${KokkosCore_C_COMPILER}")
  message(STATUS "Kokkos CMAKE_FORTRAN_COMPILER ${KokkosCore_FORTRAN_COMPILER}")

  message(STATUS "Kokkos enabled - Overriding compiler flags for consistency with Kokkos install")
  set(CMAKE_CXX_FLAGS ${KokkosCore_CXX_FLAGS})
  set(CMAKE_C_FLAGS ${KokkosCore_C_FLAGS})
  set(CMAKE_FORTRAN_FLAGS ${KokkosCore_FORTRAN_FLAGS})
  message(STATUS "Kokkos CMAKE_CXX_FLAGS ${KokkosCore_CXX_FLAGS}")
  message(STATUS "Kokkos CMAKE_C_FLAGS ${KokkosCore_C_FLAGS}")
  message(STATUS "Kokkos CMAKE_FORTRAN_FLAGS ${KokkosCore_FORTRAN_FLAGS}")
endif()

macro(engpar_export_lib target headers)
  bob_export_target(${target})
  install(FILES ${headers} DESTINATION include)
endmacro(engpar_export_lib)


option(ENABLE_PUMI "Build with PUMI interface" OFF)
if (ENABLE_PUMI)
  set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")
  find_package(SCOREC REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
  add_subdirectory(pcu)
  include_directories(pcu)
endif()


option(ENGPAR_FORTRAN_INTERFACE "enable fortran interface" OFF)
message(STATUS "ENGPAR_FORTRAN_INTERFACE: ${ENGPAR_FORTRAN_INTERFACE}")
if(ENGPAR_FORTRAN_INTERFACE)
  enable_language(Fortran)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
endif()

include_directories(${SCOREC_PREFIX}/include)


add_subdirectory(support)

add_subdirectory(agi)

#add_subdirectory(zoltan)

add_subdirectory(partition)

add_subdirectory(interfaces)

if(ENGPAR_FORTRAN_INTERFACE)
  add_subdirectory(ftn)
endif()

add_subdirectory(test)


bob_end_package()
