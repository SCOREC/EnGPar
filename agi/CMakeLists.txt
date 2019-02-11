set(PACKAGE_VERSION "1.0")


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

if (KOKKOS_ENABLED)
  target_include_directories(agi PUBLIC ${KokkosCore_INCLUDE_DIRS})
  target_link_libraries(agi PUBLIC ${KokkosCore_LIBRARIES})
endif()

target_link_libraries(agi PUBLIC engpar_support)
target_compile_definitions(agi INTERFACE ${ENGPAR_DEFINITIONS})
engpar_export_lib(agi "${AGI_HEADERS}")

bob_end_subdir()
