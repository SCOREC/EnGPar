option(ENABLE_ZOLTAN "Enable Zoltan iterface [ON|OFF]" OFF)


set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR/../cmake/}")

if (ENABLE_ZOLTAN)
  find_package(Zoltan MODULE REQUIRED)
endif()
include_directories(
  ../support
  ../agi
  ../zoltan
  ../partition
  ../interfaces
  ${ZOLTAN_INCLUDE_DIRS}
  ${PARMETIS_INCLUDE_DIRS}
  )

if(ENGPAR_FORTRAN_INTERFACE)
  set(ALL_LIBS engpar_ftn)
endif()

set(ALL_LIBS
  ${ALL_LIBS}
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
)

function(make_exe exename srcname)
  add_executable(${exename} ${srcname} buildGraphs.cpp)
  target_link_libraries(${exename} ${ALL_LIBS})
endfunction(make_exe)

make_exe(constructGraph constructGraph.cpp)
make_exe(testAdjacentTraversal testAdjacentTraversal.cpp)
make_exe(testEdgeTraversal edgeTraversal.cpp)
make_exe(testMigration testMigration.cpp)
make_exe(testTags testTags.cpp)
make_exe(splitCommunicator splitCommunicator.cpp)

make_exe(buildBinaryGraph buildBinaryGraph.cpp)

make_exe(testDistanceQueue testDistanceQueue.cpp)
make_exe(vtxBalance vtxBalance.cpp)
make_exe(balance balance.cpp)

if (ENABLE_PUMI)
  make_exe(render render.cpp)
  make_exe(testFileIO testFileIO.cpp)
  make_exe(buildMeshGraph buildMeshGraph.cpp)
  make_exe(balanceMesh balanceMesh.cpp)
endif()

if (ZOLTAN_FOUND)
  make_exe(cutVertex cutVertex.cpp)
  make_exe(pagerank pagerank.cpp)
endif()


if (ENABLE_KOKKOS)
  set(kokkos_ext "host")
  if(KOKKOS_ENABLE_CUDA)
    set(kokkos_ext "gpu")
  endif()
  make_exe(kokkosHelloWorld.${kokkos_ext} kokkosHelloWorld.cpp)
  set_property(TARGET kokkosHelloWorld.${kokkos_ext} PROPERTY CXX_STANDARD_REQUIRED 11)
  make_exe(kokkosFor.${kokkos_ext} kokkosFor.cpp)
  set_property(TARGET kokkosFor.${kokkos_ext} PROPERTY CXX_STANDARD_REQUIRED 11)
  make_exe(bfsSearch.${kokkos_ext} bfsSearch.cpp)
  set_property(TARGET bfsSearch.${kokkos_ext} PROPERTY CXX_STANDARD_REQUIRED 11)
endif()

if (ENGPAR_FORTRAN_INTERFACE)
  make_exe(ftnTest ftnTest.f90)
  set_property(TARGET ftnTest PROPERTY LINKER_LANGUAGE Fortran)
  set_source_files_properties(ftnTest.f90 PROPERTIES COMPILE_FLAGS -cpp)
endif()

#Ctest file
if(IS_TESTING)
  include(testing.cmake)
endif()
