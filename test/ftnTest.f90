program main
  use engpar
  use iso_c_binding
  implicit none
  include 'mpif.h'
#include "../agi/agi_types.h"

  integer, parameter :: NUMVERTS = 2, NUMEDGES = 1, NUMPINS = 2
  integer :: ierr, nverts, nedges, npins
  integer(AGI_GID_FT) :: gid, verts(NUMVERTS)
  integer(AGI_GID_FT) :: edges(NUMEDGES), pins(NUMPINS)
  integer(AGI_LID_FT) :: degs(NUMEDGES)
  real(AGI_WGT_FT) :: weights(NUMVERTS)
  type(c_ptr) :: graph
  logical(C_BOOL) :: isHg = .false.
  call mpi_init(ierr)
  call cengpar_initialize()
  call cengpar_setftncommunicator(MPI_COMM_WORLD)
  write (*,'(a)' ) 'Hello EnGPar!'
  write (*,* ) sizeof(gid)
  write (*, '(a,z8)' ) 'ftn pre graph', graph
  graph = cengpar_createEmptyGraph()
  write (*, '(a,z8)' ) 'ftn pre graph', graph
  nverts = NUMVERTS
  nedges = NUMEDGES
  npins = NUMPINS
  verts = (/ 0, 1 /)
  weights = (/ 1.0, 1.0 /)
  edges = (/ 0 /)
  degs = (/ 1 /)
  pins = (/ 0, 1 /)
  call cengpar_constructVerts(graph, isHg, verts, weights, nverts)
  call cengpar_constructEdges(graph, edges, degs, pins, nedges, npins)
  call cengpar_finalize()
  call mpi_finalize(ierr)
  stop
end
