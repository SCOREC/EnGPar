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

mpi_test(buildParallelBinaryGraph 4
  ./buildBinaryGraph
  "${GRAPHS}/google.ebin")

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

mpi_test(vtxBalanceCube 4
  ./vtxBalance
  "${GRAPHS}/cube/4/"
  .1)

mpi_test(vtxBalanceTorus 4
  ./vtxBalance
  "${GRAPHS}/torus/4/"
  .15)

mpi_test(vtxBalanceRing 2
  ./vtxBalance
  "${GRAPHS}/ring.ebin"
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

mpi_test(balanceTree 2
  ./balance
  "${GRAPHS}/tree.ebin")

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

if(ENGPAR_FORTRAN_INTERFACE)
  mpi_test(ftnTest 1 ./ftnTest)
endif()
