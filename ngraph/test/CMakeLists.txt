include_directories(
  ../agi
  ../interfaces
  )

set(ALL_LIBS
  ${AGI_LIBS}
  ${INTERFACE_LIBS}
  )

add_executable(buildMeshGraph buildMeshGraph.cpp)
target_link_libraries(buildMeshGraph ${ALL_LIBS})
