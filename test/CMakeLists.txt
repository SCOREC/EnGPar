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

message(${ALL_LIBS})
add_executable(buildBinaryGraph buildBinaryGraph.cpp)
target_link_libraries(buildBinaryGraph ${ALL_LIBS})
