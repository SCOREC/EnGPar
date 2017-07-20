set(PACKAGE_VERSION "1.0")
project(interfaces)


set(INTERFACE_INCLUDE_DIRS
  ../partition
  ../support
  ../agi
  )

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

include_directories(${INTERFACE_INCLUDE_DIRS} ${SCOREC_INCLUDE_DIRS})

add_library(interfaces ${INTERFACE_SOURCES})
target_link_libraries(interfaces ${DEP_LIBS})
if(ENABLE_PUMI)
  target_link_libraries(interfaces SCOREC::core)
else()
  target_link_libraries(interfaces SCOREC::pcu)
endif()

engpar_export_lib(interfaces "${INTERFACE_HEADERS}")
