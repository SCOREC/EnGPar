program main
  use engpar
  implicit none
  include 'mpif.h'

  integer :: ierr
  call MPI_INIT(ierr)
  call cengpar_initialize()
  write (*,'(a)' ) 'Hello world!'

  stop
end
