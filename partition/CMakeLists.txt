set(PACKAGE_VERSION "1.0")
project(partition)

include_directories(
  ../agi
  )

set(ENGPAR_SOURCES
  src/engpar_balancer.cpp
  src/engpar_vtxBalancer.cpp
  src/engpar_queue.cpp
  src/engpar_selector.cpp
  )

set(ENGPAR_HEADERS
  engpar.h
  )

add_library(engpar ${ENGPAR_SOURCES})

target_link_libraries(engpar PUBLIC agi)
