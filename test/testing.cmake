set(MESHES ""
    CACHE STRING
    "Extracted http://scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz")
set(GRAPHS ""
    CACHE STRING
    "Extracted http://scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz")

function(mpi_test TESTNAME PROCS EXE)
  add_test(
    NAME ${TESTNAME}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS} ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
  )
endfunction(mpi_test)

mpi_test(constructGraphSerial 1
  ./constructGraph)
mpi_test(constructGraph2 2
  ./constructGraph)
mpi_test(constructGraph4 4
  ./constructGraph)
mpi_test(testMigration2 2
  ./testMigration)
mpi_test(testMigration4 4
  ./testMigration)
mpi_test(constructAndPartition2 2
  ./constructAndPartition)
mpi_test(constructAndPartition4 4
  ./constructAndPartition)

IF(ENABLE_PUMI)
  mpi_test(buildSerialMeshGraph32 1
    ./buildMeshGraph
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi11/cube.smb
    3
    2)

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
    ${GRAPHS}/cube/cube)

  mpi_test(testFileIOParallel 4
    ./testFileIO
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi670/4/cube.smb
    ${GRAPHS}/cube/4/)

  mpi_test(testFileIOParallelTorus 4
    ./testFileIO
    ${MESHES}/torus/torus.dmg
    ${MESHES}/torus/4imb/torus.smb
    ${GRAPHS}/torus/4/
    0)

  mpi_test(testFileIOParallelTorus01 4
    ./testFileIO
    ${MESHES}/torus/torus.dmg
    ${MESHES}/torus/4imb/torus.smb
    ${GRAPHS}/torus/4_01/
    0
    1)
ENDIF()
mpi_test(buildSerialBinaryRing 1
  ./buildBinaryGraph
  "${GRAPHS}/ring.ebin")

mpi_test(buildSerialBinaryTree 1
  ./buildBinaryGraph
  "${GRAPHS}/tree.ebin")

mpi_test(buildParallelBinaryGraph2 2
  ./buildBinaryGraph
  "${GRAPHS}/gnutella.ebin")

mpi_test(buildParallelBinaryGraph4 4
  ./buildBinaryGraph
  "${GRAPHS}/enron.ebin")

mpi_test(testAdjacentSerial 1
  ./testAdjacentTraversal
  "${GRAPHS}/cube/cube"
  "${GRAPHS}/ring.ebin")

mpi_test(testAdjacentParallel 2
  ./testAdjacentTraversal
  "${GRAPHS}/cube/4/"
  "${GRAPHS}/tree.ebin")

mpi_test(testEdgeRing 1
  ./testEdgeTraversal
  "${GRAPHS}/ring.ebin")

mpi_test(testEdgeTree 2
  ./testEdgeTraversal
  "${GRAPHS}/tree.ebin")

#Diffusive Load Balancing Tests

mpi_test(testDistanceQueueSerial 1
  ./testDistanceQueue)
mpi_test(testDistanceQueueParallel 4
  ./testDistanceQueue)

# aero mesh tests
mpi_test(testDistanceQueueAero25114 1
  ./testDistanceQueue "${GRAPHS}/aero1Belm/graph128Ki_25114")
mpi_test(testDistanceQueueAero88637 1
  ./testDistanceQueue "${GRAPHS}/aero1Belm/graph128Ki_88637")
mpi_test(testDistanceQueueAero81334 1
  ./testDistanceQueue "${GRAPHS}/aero1Belm/afterMigr_81334")
mpi_test(testDistanceQueueAero1457 1
  ./testDistanceQueue "${GRAPHS}/aero1Belm/afterMigr_1457")

mpi_test(ngraphToSsgLine 1 ./ngraphToSsg)
mpi_test(ngraphToSsgAero25114 1
  ./ngraphToSsg "${GRAPHS}/aero1Belm/graph128Ki_25114")

if (ENABLE_PARMETIS)
  mpi_test(splitCube1to2 2
    ./split
    "${GRAPHS}/cube/cube"
    2
    cake/)
  
  mpi_test(splitCube1to4 4
    ./split
    "${GRAPHS}/cube/cube"
    4
    cake/)

  mpi_test(splitCube4to8 8
    ./split
    "${GRAPHS}/cube/4/"
    2
    cake/)

  mpi_test(splitTorus4to8 8
    ./split
    "${GRAPHS}/torus/4_01/"
    2
    cake/)

  mpi_test(localSplit1to4 4
    ./local_split
    "${GRAPHS}/cube/cube"
    4
    cake/)

  mpi_test(localSplit4to8 8
    ./local_split
    "${GRAPHS}/torus/4/"
    2
    cake/)

  mpi_test(splitAndBalanceCube 8
    ./splitAndBalance
    "${GRAPHS}/cube/4/"
    2
    cake/)

  mpi_test(splitAndBalanceTorus 8
    ./splitAndBalance
    "${GRAPHS}/torus/4_01/"
    2
    cake/)
  mpi_test(splitAndBalanceEqualParts 4
    ./splitAndBalance
    "${GRAPHS}/torus/4/"
    1
    cake/)

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
  
mpi_test(vtxBalanceCube 4
  ./vtxBalance
  "${GRAPHS}/cube/4/"
  .1)

mpi_test(vtxBalanceTorus 4
  ./vtxBalance
  "${GRAPHS}/torus/4/"
  .15)

mpi_test(vtxBalanceGraph 2
  ./vtxBalance
  "${GRAPHS}/gnutella.ebin"
  .1)

mpi_test(balanceCube 4
  ./balance
  "${GRAPHS}/cube/4/")

mpi_test(balanceTorus 4
  ./balance
  "${GRAPHS}/torus/4/")

mpi_test(balanceTorusVtxFaceElm 4
  ./balance
  "${GRAPHS}/torus/4_01/"
  yes)


if (ENABLE_ZOLTAN)

endif()

if (ENABLE_KOKKOS)
  mpi_test(kokkosHelloWorld 1
    ./kokkosHelloWorld)

  mpi_test(kokkosForRing 1
    ./kokkosFor
    "${GRAPHS}/ring.ebin")

  mpi_test(kokkosForTree 1
    ./kokkosFor
    "${GRAPHS}/tree.ebin")

  mpi_test(bfsSearchRing 1
    ./bfsSearch
    "${GRAPHS}/ring.ebin")

  mpi_test(bfsSearchTree 1
    ./bfsSearch
    "${GRAPHS}/tree.ebin")

endif()

if (ENGPAR_FORTRAN_INTERFACE AND ENABLE_PARMETIS)
  mpi_test(ftnTest 2 ./ftnTest)
  mpi_test(splitAndBalanceFtn 4 ./splitAndBalanceFtn
    ${GRAPHS}/cube/cube ${GRAPHS}/cube/cube_4p) 
endif()
