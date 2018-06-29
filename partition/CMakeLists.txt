include_directories(
  ./
  )

add_subdirectory(Metrics)

add_subdirectory(Diffusive)

add_subdirectory(Multilevel)

include_directories(
  Diffusive
  Multilevel
  )

set(ENGPAR_SOURCES
  engpar_input.cpp
  engpar.cpp
  )

set(ENGPAR_HEADERS
  engpar_input.h
  engpar.h
  engpar_types.h
  )

add_library(engpar ${ENGPAR_SOURCES})

target_include_directories(engpar PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(engpar PUBLIC diffusive multilevel)

engpar_export_lib(engpar "${ENGPAR_HEADERS}")

bob_end_subdir()
