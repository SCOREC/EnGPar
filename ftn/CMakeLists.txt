#TODO: Fix reading these for fortran interface
option(LONG_LOCAL_INDICES "Turn on if a 32 byte int is too small for local indices [ON|OFF]" OFF)
if (LONG_LOCAL_INDICES)
  add_definitions(-DLOCAL_ID_TYPE=int64_t)
  add_definitions(-DLOCAL_ID_FTYPE=C_INT64_T)
else()
  add_definitions(-DLOCAL_ID_TYPE=int32_t)
  add_definitions(-DLOCAL_ID_FTYPE=C_INT32_T)
endif()
option(BIG_ENDIAN "Turn on if on a big endian machine [ON|OFF]" OFF)
if (BIG_ENDIAN)
  message(STATUS "BIG_ENDIAN is enabled")
  add_definitions(-D_ENGPAR_BIG_ENDIAN)
endif()



add_library(engpar_ftn engpar.f90)
include_directories(${CMAKE_SOURCE_DIR}/agi)
set_source_files_properties(engpar.f90 PROPERTIES COMPILE_FLAGS -cpp)
set(HEADERS)
engpar_export_lib(engpar_ftn "${HEADERS}")

bob_end_subdir()
