include_directories(
  ../agi
  ../interfaces
  )

set(ALL_LIBS
  agi
  interfaces
  SCOREC::core
  )

add_executable(buildMeshGraph buildMeshGraph.cpp)
target_link_libraries(buildMeshGraph ${ALL_LIBS})
