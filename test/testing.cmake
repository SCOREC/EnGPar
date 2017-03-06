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

mpi_test(buildParallelMeshGraph31 4
  ./buildMeshGraph
  "${MESHES}/cube/cube.dmg"
  "${MESHES}/cube/pumi670/4/cube.smb"
  3
  1)

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
  "${MESHES}/cube/cube.dmg"
  "${MESHES}/cube/pumi11/cube.smb"
  "${GRAPHS}/ring.ebin")

mpi_test(testAdjacentParallel 2
  ./testAdjacentTraversal
  "${MESHES}/cube/cube.dmg"
  "${MESHES}/cube/pumi670/2/cube.smb"
  "${GRAPHS}/tree.ebin")

mpi_test(testEdgeRing 1
  ./testEdgeTraversal
  "${GRAPHS}/ring.ebin")

mpi_test(testEdgeTree 2
  ./testEdgeTraversal
  "${GRAPHS}/tree.ebin")
