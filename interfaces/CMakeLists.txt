set(PACKAGE_VERSION "1.0")
project(interfaces)


set(INTERFACE_INCLUDE_DIRS
  ../agi
  )

set(DEP_LIBS agi)

set(INTERFACE_SOURCES
  apfGraph.cpp
  binGraph.cpp
  )

set(INTERFACE_HEADERS
  apfGraph.h
  binGraph.h
  )

include_directories(${INTERFACE_INCLUDE_DIRS} ${SCOREC_INCLUDE_DIRS})

add_library(interfaces ${INTERFACE_SOURCES})
target_link_libraries(interfaces ${DEP_LIBS})
target_link_libraries(interfaces SCOREC::core)
