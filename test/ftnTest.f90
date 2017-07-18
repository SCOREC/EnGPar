program main
  use engpar
  implicit none
  include 'mpif.h'

  integer :: ierr
  call mpi_init(ierr)
  call cengpar_initialize()
  write (*,'(a)' ) 'Hello EnGPar!'
  call cengpar_finalize()
  call mpi_finalize(ierr)
  stop
end
