set(ENGPAR_SOURCES
  engpar_metrics.cpp
  )

set(ENGPAR_HEADERS
  )

add_library(engpar_metrics ${ENGPAR_SOURCES})

target_include_directories(engpar_metrics PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(engpar_metrics PUBLIC agi)

engpar_export_lib(engpar_metrics "")

bob_end_subdir()
