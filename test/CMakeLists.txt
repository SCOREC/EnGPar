include_directories(
  ../zoltan
  ${ZOLTAN_INCLUDE_DIRS}
  )

if (KOKKOS_ENABLED)
  include_directories(${KokkosCore_INCLUDE_DIRS} ${KokkosKernels_INCLUDE_DIRS})
endif()

if(ENGPAR_FORTRAN_INTERFACE)
  set(ALL_LIBS engpar_ftn)
endif()

set(ALL_LIBS
  ${ALL_LIBS}
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

if (KOKKOS_ENABLED)
  set(ALL_LIBS ${ALL_LIBS} ${KokkosCore_LIBRARIES} ${KokkosKernels_LIBRARIES})
endif()

set(TESTING_SOURCES
  buildGraphs.cpp
  gatherGraphs.cpp
  TestingSuite.cpp
  )

function(make_utility exename srcname)
  add_executable(${exename} ${srcname} buildGraphs.cpp)
  target_link_libraries(${exename} ${ALL_LIBS})
endfunction(make_utility)

function(make_test exename srcname)
  add_executable(${exename} ${srcname} ${TESTING_SOURCES})
  target_link_libraries(${exename} ${ALL_LIBS})
endfunction(make_test)

make_utility(addCoords addCoords.cpp)
make_utility(balance balance.cpp)
make_utility(split split.cpp)
make_utility(splitAndBalance splitAndBalance.cpp)

if (ENABLE_PUMI)
  make_utility(render render.cpp)
  make_utility(testFileIO testFileIO.cpp)
  make_utility(buildMeshGraph buildMeshGraph.cpp)
  make_utility(balanceMesh balanceMesh.cpp)
  make_utility(splitAndBalanceMesh splitAndBalanceMesh.cpp)
endif()

if (ZOLTAN_FOUND)
  make_utility(cutVertex cutVertex.cpp)
  make_utility(pagerank pagerank.cpp)
endif()


if (KOKKOS_ENABLED)
  make_utility(kokkosHelloWorld kokkosHelloWorld.cpp)
  make_utility(bfsSearch bfsSearch.cpp)
  make_utility(kokkosColoring kokkosColoring.cpp)
endif()

if (ENGPAR_FORTRAN_INTERFACE AND ENABLE_PARMETIS)
  make_utility(ftnTest ftnTest.f90)
  set_property(TARGET ftnTest PROPERTY LINKER_LANGUAGE Fortran)
  set_source_files_properties(ftnTest.f90 PROPERTIES COMPILE_FLAGS -cpp)
  make_utility(splitAndBalanceFtn splitAndBalance.f90)
  set_property(TARGET splitAndBalanceFtn PROPERTY LINKER_LANGUAGE Fortran)
  set_source_files_properties(splitAndBalance.f90 PROPERTIES COMPILE_FLAGS -cpp)
  make_utility(splitFtn split.f90)
  set_property(TARGET splitFtn PROPERTY LINKER_LANGUAGE Fortran)
  set_source_files_properties(split.f90 PROPERTIES COMPILE_FLAGS -cpp)
endif()

#Ctest file
if(IS_TESTING)
  set(MESHES ""
    CACHE STRING
    "See sub module pumi_meshes")
  set(GRAPHS ""
    CACHE STRING
    "See sub module EnGPar-graphs")
  add_definitions(-DENGPAR_GRAPHS="${GRAPHS}")
  make_test(ConstructTestSuite ConstructTestSuite.cpp)
  make_test(NgraphTestSuite NgraphTestSuite.cpp)
  make_test(PartitionTestSuite PartitionTestSuite.cpp)
  include(testing.cmake)
endif()

bob_end_subdir()
