set(ENGPAR_SOURCES
  engpar_diffusive_input.cpp
  engpar_weight_input.cpp
  src/engpar_balancer.cpp
  src/engpar_vtxBalancer.cpp
  src/engpar_weightBalancer.cpp
  src/engpar_queue.cpp
  src/engpar_selector.cpp
  src/engpar_weightSelector.cpp
  src/engpar_make.cpp
  src/engpar_sd.cpp
  )

set(ENGPAR_HEADERS
  engpar_diffusive_input.h
  engpar_weight_input.h
  )

add_library(diffusive ${ENGPAR_SOURCES})

target_include_directories(diffusive PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(diffusive PUBLIC engpar_metrics)

engpar_export_lib(diffusive "${ENGPAR_HEADERS}")

bob_end_subdir()
