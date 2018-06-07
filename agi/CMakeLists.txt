set(PACKAGE_VERSION "1.0")
project(agi)

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
  )

add_library(agi ${AGI_SOURCES})

target_include_directories(agi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

if (ENABLE_KOKKOS)
  target_link_libraries(agi PUBLIC kokkos)
endif()
target_link_libraries(agi PUBLIC engpar_support)

engpar_export_lib(agi "${AGI_HEADERS}")
