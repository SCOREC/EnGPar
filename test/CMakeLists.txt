include_directories(
  ../agi
  ../zoltan
  ../interfaces
  )
set(ALL_LIBS
  agi
  interfaces
  agi_zoltan
  SCOREC::core
  )

add_executable(buildMeshGraph buildMeshGraph.cpp)
target_link_libraries(buildMeshGraph ${ALL_LIBS})


add_executable(buildBinaryGraph buildBinaryGraph.cpp)
target_link_libraries(buildBinaryGraph ${ALL_LIBS})

option(ENABLE_ZOLTAN "Enable Zoltan itergace [ON|OFF]" OFF)
message(STATUS "ENABLE_ZOLTAN: " ${ENABLE_ZOLTAN})

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR/../cmake/}")

if (ENABLE_ZOLTAN)
  find_package(Zoltan MODULE REQUIRED)
  add_executable(cutVertex cutVertex.cpp)
  target_link_libraries(cutVertex PUBLIC ${ALL_LIBS})

  target_include_directories(cutVertex PRIVATE
    ${ZOLTAN_INCLUDE_DIRS}
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(cutVertex PUBLIC
    ${ZOLTAN_LIBRARIES}
    ${PARMETIS_LIBRARIES}
    )
  add_executable(pagerank pagerank.cpp)
  target_link_libraries(pagerank PUBLIC ${ALL_LIBS})

  target_include_directories(pagerank PRIVATE
    ${ZOLTAN_INCLUDE_DIRS}
    ${PARMETIS_INCLUDE_DIRS}
    )
  target_link_libraries(pagerank PUBLIC
    ${ZOLTAN_LIBRARIES}
    ${PARMETIS_LIBRARIES}
    )
endif()

if (ENABLE_KOKKOS)
  find_package(Kokkos MODULE REQUIRED)
  add_executable(kokkosHelloWorld kokkosHelloWorld.cpp)
  target_link_libraries(kokkosHelloWorld PUBLIC ${ALL_LIBS})
  target_include_directories(kokkosHelloWorld PRIVATE ${KOKKOS_INCLUDE_DIRS})
  target_link_libraries(kokkosHelloWorld PUBLIC ${KOKKOS_LIBRARIES})
endif()

#Ctest file
if(IS_TESTING)
  include(testing.cmake)
endif()
