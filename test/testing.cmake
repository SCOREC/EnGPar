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
ENDIF()

if (ENABLE_PARMETIS)
  mpi_test(splitCube1to2 2
    ./split
    "${GRAPHS}/cube/cube"
    1
    cake/)
  
  mpi_test(splitCube1to4 4
    ./split
    "${GRAPHS}/cube/cube"
    1
    cake/)

  mpi_test(splitCube4to8 8
    ./split
    "${GRAPHS}/cube/4/"
    4
    cake/)

  mpi_test(splitTorus4to6 6
    ./split
    "${GRAPHS}/torus/4/"
    4
    cake/)

  mpi_test(splitTorus4to8 8
    ./split
    "${GRAPHS}/torus/4_01/"
    4
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

  mpi_test(balanceMultipleTorus 4
    ./balanceMultiple
    "${GRAPHS}/torus/4/")

  mpi_test(balanceWeights4_100 4
    ./balanceWeights
    100
    1.05)

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
    ${GRAPHS}/cube/cube ${GRAPHS}/cube/cube_4p 4)
  mpi_test(splitFtn 5 ./splitFtn
    ${GRAPHS}/torus/4/ ${GRAPHS}/torus/5/ 4)
endif()
