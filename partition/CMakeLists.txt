include_directories(
  ./
  )

add_subdirectory(Metrics)

add_subdirectory(Diffusive)

add_subdirectory(Multilevel)

add_subdirectory(Coloring)

include_directories(
  Diffusive
  Multilevel
  Coloring
  )

set(ENGPAR_SOURCES
  engpar_input.cpp
  engpar.cpp
  )

set(ENGPAR_HEADERS
  engpar_input.h
  engpar.h
  engpar_ptn_types.h
  )


add_library(engpar ${ENGPAR_SOURCES})

target_include_directories(engpar PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(engpar PUBLIC agi diffusive multilevel coloring)

if (ENABLE_PARMETIS)
  target_include_directories(engpar PRIVATE
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(engpar PUBLIC
    ${PARMETIS_LIBRARIES}
    )
endif()

engpar_export_lib(engpar "${ENGPAR_HEADERS}")

bob_end_subdir()
