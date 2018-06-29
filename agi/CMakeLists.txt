set(PACKAGE_VERSION "1.0")

set(AGI_DEFINITIONS
  )
option(LONG_LOCAL_INDICES "Turn on if a 32 byte int is too small for local indices [ON|OFF]" OFF)
if (LONG_LOCAL_INDICES)
  message(STATUS "Local indices are 64 bytes")
  set(AGI_DEFINITIONS
    ${AGI_DEFINITIONS}
    -DLOCAL_ID_TYPE=int64_t
    -DLOCAL_ID_FTYPE=C_INT64_T)
  add_definitions(-DLOCAL_ID_TYPE=int64_t)
  add_definitions(-DLOCAL_ID_FTYPE=C_INT64_T)
else()
  message(STATUS "Local indices are 32 bytes")
    set(AGI_DEFINITIONS
      ${AGI_DEFINITIONS}
      -DLOCAL_ID_TYPE=int32_t
      -DLOCAL_ID_FTYPE=C_INT32_T)
  add_definitions(-DLOCAL_ID_TYPE=int32_t)
  add_definitions(-DLOCAL_ID_FTYPE=C_INT32_T)
endif()
option(BIG_ENDIAN "Turn on if on a big endian machine [ON|OFF]" OFF)
if (BIG_ENDIAN)
  message(STATUS "BIG_ENDIAN is enabled")
  set(AGI_DEFINITIONS
    ${AGI_DEFINITIONS}
    -DENGPAR_BIG_ENDIAN)
  add_definitions(-D_ENGPAR_BIG_ENDIAN)
endif()


set(AGI_SOURCES
  agiConstruct.cpp
  agiVertex.cpp
  agiEdge.cpp
  agiIterate.cpp
  agiTag.cpp
  agiMigration.cpp
  agiWeightMigration.cpp
  agiMigrationTimers.cpp
  agiMigrate.cpp
  agiWeightMigrate.cpp
  agiIO.cpp
  agiValid.cpp
  agiVTK.cpp
  agi_typeconvert.cpp
  )

set(AGI_HEADERS
  ngraph.h
  pngraph.h
  agi.h
  agi_types.h
  
  )

add_library(agi ${AGI_SOURCES})

target_include_directories(agi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

if (ENABLE_KOKKOS)
  set(AGI_DEFINITIONS
    ${AGI_DEFINITIONS}
    -DKOKKOS_ENABLED)
  add_definitions(-DKOKKOS_ENABLED)
  target_link_libraries(agi PUBLIC kokkos)
endif()
target_link_libraries(agi PUBLIC engpar_support)
target_compile_definitions(agi INTERFACE ${AGI_DEFINITIONS})
engpar_export_lib(agi "${AGI_HEADERS}")

bob_end_subdir()
