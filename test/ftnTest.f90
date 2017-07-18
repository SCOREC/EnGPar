program main
  use engpar
  use iso_c_binding
  implicit none
  include 'mpif.h'
#include "../agi/agi_types.h"

  integer :: ierr
  integer(AGI_GID_FT) :: gid
  type(c_ptr) :: graph
  call mpi_init(ierr)
  call cengpar_initialize()
  call cengpar_setftncommunicator(MPI_COMM_WORLD)
  write (*,'(a)' ) 'Hello EnGPar!'
  write (*,* ) sizeof(gid)
  call cengpar_finalize()
  call mpi_finalize(ierr)
  stop
end
