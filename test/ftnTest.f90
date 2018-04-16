program main
  use engpar
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer, parameter :: nghosts = 1
  integer :: nverts, nedges, npins
  integer :: ierr, self
  integer(ENGPAR_GID_T), dimension(:), allocatable :: verts, edges, pins
  integer(ENGPAR_LID_T), dimension(:), allocatable :: degs
  real(ENGPAR_WGT_T), dimension(:), allocatable :: weights, eweights
  integer(ENGPAR_GID_T), dimension(nghosts) :: ghostverts
  integer(ENGPAR_PART_T), dimension(nghosts):: ghostowners
  integer(ENGPAR_PART_T), dimension(:), allocatable :: parts
  type(c_ptr) :: graph
  logical(C_BOOL) :: isHg = .false.
  real(C_DOUBLE) :: tol, stepfactor
  integer :: verbosity
  integer :: i
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call cengpar_initialize()
  call cengpar_setftncommunicator(MPI_COMM_WORLD)
  graph = cengpar_createEmptyGraph()
  !  p0   |   p1
  ! 1-2-3-|-4-5-6-7-8
  if ( self == 0 ) then
    ! 1-2-3
    nverts = 3
    nedges = 3
    npins = nedges*2
  else
    ! 4-5-6-7-8
    nverts = 5
    nedges = 5
    npins = nedges*2
  end if
  allocate(verts(nverts))
  allocate(weights(nverts))
  allocate(edges(nedges))
  allocate(pins(npins))
  allocate(degs(nedges))
  allocate(eweights(nedges))
  if ( self == 0 ) then
    ! constructVerts args
    verts   = (/ 1,2,3 /)
    weights = (/ 1.0,1.0,1.0 /)
    ! constructEdges args
    edges   = (/ 1,2,3 /)
    degs    = (/ 2,2,2 /)
    eweights = (/ 1.0,1.0,1.0 /)
    pins    = (/ 1,2, 2,3, 3,4 /)
    ! construcGhosts args
    ghostverts = (/ 4 /)
    ghostowners = (/ 1 /)
  else
    ! constructVerts args
    verts   = (/ 4,5,6,7,8 /)
    weights = (/ 1.0,1.0,1.0,1.0,1.0 /)
    ! constructEdges args
    edges   = (/ 3,4,5,6,7 /)
    degs    = (/ 2,2,2,2,2 /)
    eweights = (/ 1.0,1.0,1.0,1.0,1.0 /)
    pins    = (/ 4,3, 4,5, 5,6, 6,7, 7,8 /)
    ! construcGhosts args
    ghostverts = (/ 3 /)
    ghostowners = (/ 0 /)
  end if
  call cengpar_constructVerts(graph, isHg, verts, weights, nverts)
  call cengpar_constructEdges(graph, edges, degs, eweights, pins, nedges, npins)
  call cengpar_constructGhosts(graph, ghostverts, ghostowners, nghosts)
  call cengpar_checkValidity(graph);
  tol = 1.05
  stepfactor = 0.1
  verbosity = 1
  call cengpar_balanceVertices(graph, tol, stepfactor, verbosity);
  allocate(parts(nverts))
  ! verts and parts are output arrays, verts[i] is assigned to process
  !   parts[i], i=1:nverts
  call cengpar_getPartition(graph, verts, parts, nverts)
  write(*,*) '<idx>, <vertex id>, <process id>'
  do i=1, nverts
    write(*,*) i, verts(i), parts(i)
  end do
  call cengpar_destroyGraph(graph);
  call cengpar_finalize()
  call mpi_finalize(ierr)
  stop
end
