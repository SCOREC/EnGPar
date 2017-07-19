program main
  use engpar
  use iso_c_binding
  implicit none
  include 'mpif.h'
#include "../agi/agi_types.h"

  integer, parameter :: NUMVERTS = 2
  integer :: ierr, nverts
  integer(AGI_GID_FT) :: gid, verts(NUMVERTS)
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
  verts = (/ 0, 1 /)
  weights = (/ 1.0, 1.0 /)
  call cengpar_constructVerts(graph, isHg, verts, weights, nverts)
  call cengpar_finalize()
  call mpi_finalize(ierr)
  stop
end
