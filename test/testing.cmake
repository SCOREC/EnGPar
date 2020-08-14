set(MESHES ""
    CACHE STRING
    "See sub module pumi_meshes")
set(GRAPHS ""
    CACHE STRING
    "See sub module EnGPar-graphs")

function(mpi_test TESTNAME PROCS EXE)
  add_test(
    NAME ${TESTNAME}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS} ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
  )
endfunction(mpi_test)

mpi_test(ConstructTest_1 1
  ./ConstructTestSuite)
mpi_test(ConstructTest_2 2
  ./ConstructTestSuite)
mpi_test(ConstructTest_4 4
  ./ConstructTestSuite)

mpi_test(NgraphTest_1 1
  ./NgraphTestSuite)
mpi_test(NgraphTest_2 2
  ./NgraphTestSuite)
mpi_test(NgraphTest_4 4
  ./NgraphTestSuite)


IF(ENABLE_PUMI)
  mpi_test(buildSerialMeshGraph32 1
    ./buildMeshGraph
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi11/cube.smb
    3
    2)

  mpi_test(buildSerialMeshGraph01 1
    ./buildMeshGraph
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi11/cube.smb
    0
    1)

  mpi_test(buildSerialMeshGraph13 1
    ./buildMeshGraph
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi11/cube.smb
    1
    3)

  mpi_test(buildSerialMeshGraph31 1
    ./buildMeshGraph
    "${MESHES}/cube/cube.dmg"
    "${MESHES}/cube/pumi11/cube.smb"
    3
    1)

  mpi_test(buildParallelMeshGraph32 2
    ./buildMeshGraph
    "${MESHES}/cube/cube.dmg"
    "${MESHES}/cube/pumi670/2/cube.smb"
    3
    2)

  mpi_test(buildParallelMeshGraph30 4
    ./buildMeshGraph
    "${MESHES}/torus/torus.dmg"
    "${MESHES}/torus/4imb/torus.smb"
    3
    0)

  mpi_test(testFileIOSerial 1
    ./testFileIO
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi670/cube.smb
    ${GRAPHS}/cube/cube
    3)

  mpi_test(testFileIOParallel 4
    ./testFileIO
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi670/4/cube.smb
    ${GRAPHS}/cube/4/
    3)

  mpi_test(testFileIOParallelTorus 4
    ./testFileIO
    ${MESHES}/torus/torus.dmg
    ${MESHES}/torus/4imb/torus.smb
    ${GRAPHS}/torus/4/
    3
    0)

  mpi_test(testFileIOParallelTorus01 4
    ./testFileIO
    ${MESHES}/torus/torus.dmg
    ${MESHES}/torus/4imb/torus.smb
    ${GRAPHS}/torus/4_01/
    3
    0
    1)
  
  mpi_test(renderCube 4
    ./render
    "${GRAPHS}/cube/4/"
    )
  mpi_test(renderTorus 4
    ./render
    "${GRAPHS}/torus/4/"
    )
ENDIF()

mpi_test(VertexBalancer_2 2
  ./PartitionTestSuite 0)
mpi_test(VertexBalancer_4 4
  ./PartitionTestSuite 0)
mpi_test(MultiCriteriaBalancer_2 2
  ./PartitionTestSuite 1)
mpi_test(MultiCriteriaBalancer_4 4
  ./PartitionTestSuite 1)
mpi_test(MultipleBalances_2 2
  ./PartitionTestSuite 2)
mpi_test(MultipleBalances_4 4
  ./PartitionTestSuite 2)
mpi_test(WeightBalancer_2 2
  ./PartitionTestSuite 3)
mpi_test(WeightBalancer_4 4
  ./PartitionTestSuite 3)
mpi_test(PartWeightBalancer_2 2
  ./PartitionTestSuite 4)
mpi_test(PartWeightBalancer_4 4
  ./PartitionTestSuite 4)
if (ENABLE_PARMETIS)
  mpi_test(GlobalSplit_2 2
    ./PartitionTestSuite 10)
  mpi_test(GlobalSplit_4 4
    ./PartitionTestSuite 10)
  mpi_test(LocalSplit_2 2
    ./PartitionTestSuite 11)
  mpi_test(LocalSplit_4 4
    ./PartitionTestSuite 11)
  mpi_test(SplitAndBalance_2 2
    ./PartitionTestSuite 13)
  mpi_test(SplitAndBalance_4 4
    ./PartitionTestSuite 13)
  mpi_test(SplitAndBalance_8 8
    ./PartitionTestSuite 13)

  IF(ENABLE_PUMI)
    mpi_test(splitAndBalanceMeshPUMI 8
      ./splitAndBalanceMesh
      "${MESHES}/torus/torus.dmg"
      "${MESHES}/torus/4imb/torus.smb"
      0
      2
      )
    mpi_test(splitAndBalanceMeshEnGPar 8
      ./splitAndBalanceMesh
      "${MESHES}/torus/torus.dmg"
      "${MESHES}/torus/4imb/torus.smb"
      1
      2
      )
  endif()
endif()

if (ENABLE_ZOLTAN)
  mpi_test(SplitPHG_2 2
    ./PartitionTestSuite 12)
  mpi_test(SplitPHG_4 4
    ./PartitionTestSuite 12)
  mpi_test(SplitPHG_8 8
    ./PartitionTestSuite 12)
endif()

if (ENABLE_KOKKOS)
  mpi_test(kokkosHelloWorld 1
    ./kokkosHelloWorld)

  mpi_test(bfsSearchRing 1
    ./bfsSearch
    "${GRAPHS}/ring.ebin")

  mpi_test(bfsSearchTree 1
    ./bfsSearch
    "${GRAPHS}/tree.ebin")

  mpi_test(colorRing 1
    ./kokkosColoring
    "${GRAPHS}/ring.ebin")

  mpi_test(colorTree 1
    ./kokkosColoring
    "${GRAPHS}/tree.ebin")

  mpi_test(colorEnron 1
    ./kokkosColoring
    "${GRAPHS}/enron.ebin")

  mpi_test(colorNutella 1
    ./kokkosColoring
    "${GRAPHS}/gnutella.ebin")

  mpi_test(colorAsic 1
    ./kokkosColoring
    "${GRAPHS}/asic.ebin")

  mpi_test(colorMeshPipe 1
    ./kokkosMeshColoring
    "${MESHES}/pipe/pipe.dmg"
    "${MESHES}/pipe/pipe.smb"
    1)

  mpi_test(colorMeshTorus 4
    ./kokkosMeshColoring
    "${MESHES}/torus/torus.dmg"
    "${MESHES}/torus/4imb/torus.smb"
    1)
  
  mpi_test(colorMeshUpright 1
    ./kokkosMeshColoring
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/67k.smb"
    1)

  mpi_test(bfsMeshSquare_p2 2
    ./kokkosMeshBfs
    "${MESHES}/square/square.dmg"
    "${MESHES}/square/2/")

  mpi_test(bfsMeshMpasOcean_p2 2
    ./kokkosMeshBfs
    "${MESHES}/mpas/ocean.dmg"
    "${MESHES}/mpas/2ref/")

  mpi_test(bfsMeshUpright 4
    ./kokkosMeshBfs
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/67k_4/")

  mpi_test(bfsMeshUprightBig2 2
    ./kokkosMeshBfs
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/1.6M/2p/")

  mpi_test(bfsMeshUprightBig4 4
    ./kokkosMeshBfs
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/1.6M/4p/")

  mpi_test(bfsMeshUpright13M_2 2
    ./kokkosMeshBfs
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/13M/2p/")

  mpi_test(bfsMeshUpright13M_4 4
    ./kokkosMeshBfs
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/13M/4p/")

  mpi_test(bfsMeshUpright13M_8 8
    ./kokkosMeshBfs
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/13M/8p/")

  mpi_test(splitAndBalanceMeshEnGPar_1M_2 2
    ./splitAndBalanceMesh
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/1.6M/2p/"
    0
    1
    )

  mpi_test(splitAndBalanceMeshEnGPar_1M_4 4
    ./splitAndBalanceMesh
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/1.6M/4p/"
    0
    1
    )

  mpi_test(splitAndBalanceMeshEnGPar_13M_2 2
    ./splitAndBalanceMesh
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/13M/2p/"
    0
    1
    )

  mpi_test(splitAndBalanceMeshEnGPar_13M_4 4
    ./splitAndBalanceMesh
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/13M/4p/"
    0
    1
    )

  mpi_test(splitAndBalanceMeshEnGPar_13M_8 8
    ./splitAndBalanceMesh
    "${MESHES}/upright/upright.dmg"
    "${MESHES}/upright/13M/8p/"
    0
    1
    )

endif()

if (ENGPAR_FORTRAN_INTERFACE AND ENABLE_PARMETIS)
  mpi_test(ftnTest 2 ./ftnTest)
  mpi_test(splitAndBalanceFtn 4 ./splitAndBalanceFtn
    ${GRAPHS}/cube/cube ${GRAPHS}/cube/cube_4p 4)
  mpi_test(splitFtn 5 ./splitFtn
    ${GRAPHS}/torus/4/ ${GRAPHS}/torus/5/ 4)
endif()
