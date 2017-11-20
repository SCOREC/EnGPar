set(PACKAGE_VERSION "1.0")
project(interfaces)

set(DEP_LIBS agi)

set(INTERFACE_SOURCES
  binGraph.cpp
  cengpar.cpp
  )

set(INTERFACE_HEADERS
  binGraph.h
  cengpar.h
  )

if(ENABLE_PUMI)
  set(INTERFACE_HEADERS apfGraph.h ${INTERFACE_HEADERS})
  set(INTERFACE_SOURCES apfGraph.cpp ${INTERFACE_SOURCES})
endif()

add_library(interfaces ${INTERFACE_SOURCES})

target_include_directories(agi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(interfaces ${DEP_LIBS})

if(ENABLE_PUMI)
  target_link_libraries(interfaces SCOREC::core)
else()
  target_link_libraries(interfaces SCOREC::pcu)
endif()

engpar_export_lib(interfaces "${INTERFACE_HEADERS}")
