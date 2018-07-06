set(COLORING_SOURCE engpar_kokkosColoring.cpp)

add_library(coloring ${COLORING_SOURCE})

target_link_libraries(coloring PUBLIC agi)

target_include_directories(coloring PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

if (ENABLE_KOKKOS)
  target_include_directories(coloring PUBLIC
    ${KokkosCore_INCLUDE_DIRS} ${KokkosKernels_INCLUDE_DIRS}
    )
  target_link_libraries(coloring PUBLIC
    ${KokkosCore_LIBRARIES} ${KokkosKernels_LIBRARIES}
    )
endif()

engpar_export_lib(coloring "")
