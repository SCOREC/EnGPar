set(COLORING_SOURCE engpar_kokkosColoring.cpp)
set(COLORING_HEADERS engpar_coloring_input.h)

add_library(coloring ${COLORING_SOURCE})

target_link_libraries(coloring PUBLIC agi)

target_include_directories(coloring PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

if (ENABLE_KOKKOS)
  target_link_libraries(coloring PUBLIC Kokkos::kokkos Kokkos::kokkoskernels)
endif()

engpar_export_lib(coloring ${COLORING_HEADERS})

bob_end_subdir()
