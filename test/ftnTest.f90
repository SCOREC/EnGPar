program main
  use engpar
  implicit none
  include 'mpif.h'

  integer :: ierr
  call mpi_init(ierr)
  call cengpar_initialize()
  call cengpar_setftncommunicator(MPI_COMM_WORLD)
  write (*,'(a)' ) 'Hello EnGPar!'
  call cengpar_finalize()
  call mpi_finalize(ierr)
  stop
end
