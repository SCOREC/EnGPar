set(PACKAGE_VERSION "1.0")
project(partition)

include_directories(
  ../support
  ../agi
  )

set(ENGPAR_SOURCES
  engpar.cpp
  src/engpar_balancer.cpp
  src/engpar_vtxBalancer.cpp
  src/engpar_queue.cpp
  src/engpar_selector.cpp
  src/engpar_input.cpp
  src/engpar_make.cpp
  src/engpar_sd.cpp
  )

set(ENGPAR_HEADERS
  engpar.h
  )

add_library(engpar ${ENGPAR_SOURCES})

target_link_libraries(engpar PUBLIC agi)

engpar_export_lib(engpar "${ENGPAR_HEADERS}")
