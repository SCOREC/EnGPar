option(ENABLE_ZOLTAN "Enable Zoltan iterface [ON|OFF]" OFF)


set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR/../cmake/}")

if (ENABLE_ZOLTAN)
  find_package(Zoltan MODULE REQUIRED)
endif()
if (ENABLE_KOKKOS)
  find_package(Kokkos MODULE REQUIRED)
endif()
include_directories(
  ../support
  ../agi
  ../zoltan
  ../partition
  ../interfaces
  ${ZOLTAN_INCLUDE_DIRS}
  ${PARMETIS_INCLUDE_DIRS}
  ${KOKKOS_INCLUDE_DIRS}
  )

set(ALL_LIBS
  engpar_support
  agi
  interfaces
  agi_zoltan
  engpar
)

if(ENABLE_PUMI)
  set(ALL_LIBS ${ALL_LIBS} SCOREC::core)
else()
  set(ALL_LIBS ${ALL_LIBS} SCOREC::pcu)
endif()

set(ALL_LIBS ${ALL_LIBS}
  ${ZOLTAN_LIBRARIES}
  ${PARMETIS_LIBRARIES}
  ${KOKKOS_LIBRARIES}
)

function(make_exe exename srcname)
  add_executable(${exename} ${srcname} buildGraphs.cpp)
  target_link_libraries(${exename} ${ALL_LIBS})
endfunction(make_exe)

make_exe(constructGraph constructGraph.cpp)
make_exe(testMigration testMigration.cpp)
make_exe(testTags testTags.cpp)
make_exe(splitCommunicator splitCommunicator.cpp)
make_exe(buildBinaryGraph buildBinaryGraph.cpp)
make_exe(testEdgeTraversal edgeTraversal.cpp)
if (ENABLE_PUMI)
  make_exe(testFileIO testFileIO.cpp)
  make_exe(vtxBalance vtxBalance.cpp)
  make_exe(buildMeshGraph buildMeshGraph.cpp)
  make_exe(balance balance.cpp)
  make_exe(testAdjacentTraversal testAdjacentTraversal.cpp)
endif()

if (ZOLTAN_FOUND)
  make_exe(cutVertex cutVertex.cpp)
  make_exe(pagerank pagerank.cpp)
endif()


if (KOKKOS_FOUND)
  make_exe(kokkosHelloWorld kokkosHelloWorld.cpp)
  make_exe(kokkosFor kokkosFor.cpp)
  make_exe(bfsSearch bfsSearch.cpp)
endif()

#Ctest file
if(IS_TESTING)
  include(testing.cmake)
endif()
