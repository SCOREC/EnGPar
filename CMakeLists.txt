#Setting Version Number, Project Name
cmake_minimum_required (VERSION 3.0)
project(EnGPar VERSION 1.1.0 LANGUAGES CXX C)

include(cmake/bob.cmake)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

bob_begin_package()

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_DIRS})

find_package(Torch REQUIRED)

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
set(MPIRUN_PROCFLAG "-np" CACHE STRING
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")
set(VALGRIND "" CACHE STRING
  "the valgrind executable")
set(VALGRIND_ARGS "" CACHE STRING
  "the command line arguments to VALGRIND")

option(ENGPAR_USE_CUDA "Build for CUDA" OFF)
option(ENGPAR_USE_OPENMP "Build for OpenMP" OFF)


if (ENABLE_KOKKOS)
  find_package(Kokkos REQUIRED)
  find_package(KokkosKernels REQUIRED)
  add_definitions(-DKOKKOS_ENABLED)
  set(KOKKOS_ENABLED true)
  include_directories(${KokkosCore_INCLUDE_DIRS})

  message(STATUS "Kokkos enabled - Overriding compiler selection for consistency with Kokkos install")
  set(CMAKE_CXX_COMPILER ${KokkosCore_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${KokkosCore_C_COMPILER})
  message(STATUS "Kokkos CMAKE_CXX_COMPILER ${KokkosCore_CXX_COMPILER}")
  message(STATUS "Kokkos CMAKE_C_COMPILER ${KokkosCore_C_COMPILER}")

  message(STATUS "Kokkos enabled - Overriding compiler flags for consistency with Kokkos install")
  set(CMAKE_CXX_FLAGS ${KokkosCore_CXX_FLAGS})
  set(CMAKE_C_FLAGS ${KokkosCore_C_FLAGS})
  message(STATUS "Kokkos CMAKE_CXX_FLAGS ${KokkosCore_CXX_FLAGS}")
  message(STATUS "Kokkos CMAKE_C_FLAGS ${KokkosCore_C_FLAGS}")

  include(${CMAKE_CURRENT_LIST_DIR}/cmake/detect_trilinos_opts.cmake)
  detect_trilinos_opts()
  if(KokkosCore_HAS_CUDA AND (NOT KokkosCore_HAS_CUDA_LAMBDA))
    message(FATAL_ERROR
      "Please reconfigure Kokkos with
      -DKokkos_ENABLE_Cuda_Lambda:BOOL=ON")
  endif()
  if(ENGPAR_USE_CUDA AND (NOT KokkosCore_HAS_CUDA))
    message(FATAL_ERROR
      "Kokkos was not build with CUDA, reconfigure EnGPar"
      "-DENGPAR_USE_CUDA:BOOL=OFF")
  endif()
endif()

set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")

option(ENABLE_PUMI "Build with PUMI interface" OFF)

if (ENABLE_PUMI)
  find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
  find_package(SCOREC 2.1.0 REQUIRED COMPONENTS pcu CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
endif()


option(ENGPAR_FORTRAN_INTERFACE "enable fortran interface" OFF)
message(STATUS "ENGPAR_FORTRAN_INTERFACE: ${ENGPAR_FORTRAN_INTERFACE}")
if(ENGPAR_FORTRAN_INTERFACE)
  enable_language(Fortran)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
endif()

include_directories(${SCOREC_PREFIX}/include)

macro(engpar_export_lib target headers)
  bob_export_target(${target})
  install(FILES ${headers} DESTINATION include)
endmacro(engpar_export_lib)

set(ENGPAR_DEFINITIONS)
if(ENGPAR_USE_CUDA)
  message(STATUS "EnGPar using CUDA")
  set(ENGPAR_DEFINITIONS ${ENGPAR_DEFINITIONS} -DENGPAR_USE_CUDA)
endif()
if(ENGPAR_USE_OPENMP)
  message(STATUS "EnGPar using OpenMP")
  set(ENGPAR_DEFINITIONS ${ENGPAR_DEFINITIONS} -DENGPAR_USE_OPENMP)
endif()
if(ENGPAR_USE_OPENMP AND ENGPAR_USE_CUDA)
  message(FATAL_ERROR "EnGPar cannot use"
    "OpenMP and CUDA at the same time")
endif()

option(LONG_LOCAL_INDICES "Turn on if a 32 byte int is too small for local indices [ON|OFF]" OFF)
if (LONG_LOCAL_INDICES)
  message(STATUS "Local indices are 64 bytes")
  set(ENGPAR_DEFINITIONS
    ${ENGPAR_DEFINITIONS}
    -DLOCAL_ID_TYPE=int64_t
    -DLOCAL_ID_FTYPE=C_INT64_T)
else()
  message(STATUS "Local indices are 32 bytes")
    set(ENGPAR_DEFINITIONS
      ${ENGPAR_DEFINITIONS}
      -DLOCAL_ID_TYPE=int32_t
      -DLOCAL_ID_FTYPE=C_INT32_T)
endif()
option(BIG_ENDIAN "Turn on if on a big endian machine [ON|OFF]" OFF)
if (BIG_ENDIAN)
  message(STATUS "BIG_ENDIAN is enabled")
  set(ENGPAR_DEFINITIONS
    ${ENGPAR_DEFINITIONS}
    -DENGPAR_BIG_ENDIAN)
endif(BIG_ENDIAN)

foreach(def ${ENGPAR_DEFINITIONS})
  add_definitions(${def})
endforeach()

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
