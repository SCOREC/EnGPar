subroutine switchToOriginals(smallSize, inSmall, newComm)
  use engpar
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: smallSize, verbosity
  logical(C_BOOL) :: inSmall
  integer(kind=4) :: newComm
  integer :: ierr, self, group, groupRank
  integer :: mpiErr, strLen
  character(len=1024) :: errStr
  
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)

  if (self .lt. smallSize) then
    inSmall = .true.
  end if

  if (inSmall) then
    group=0
    groupRank=self
  else
    group = 1
    groupRank = 0
  end if
  call mpi_comm_split(MPI_COMM_WORLD, group, groupRank, newComm, ierr)
  if (ierr .ne. MPI_SUCCESS) then
    call mpi_error_string(ierr, errStr, strLen, mpiErr)
    write (*,*) "ERROR mpi comm split failed with error ierr", errStr, "!... exiting"
    stop
  end if
end subroutine

module splitHelperFns
  contains
  function cleanString(str)
    use iso_c_binding
    implicit none
    character(len=256) :: cleanString
    character(len=256) :: str
    cleanString = trim(str)//c_null_char
  end function
  subroutine getArgs(self, inGraph, outGraph, smallSize)
    implicit none
    integer :: self, numArgs, strlen, smallSize
    character(len=256) :: inGraph, outGraph, arg
    numArgs = command_argument_count()
    if ( numArgs .ne. 3 ) then
      if (self==0) write(*,*) "Usage: splitNftn </path/to/input/graph> </path/to/output/graph> <inputNumberOfParts>"
      stop
    end if

    call get_command_argument(1,inGraph)
    call get_command_argument(2,outGraph)
    call get_command_argument(3,arg)
    read(arg,*) smallSize

    if (self==0) then
      write(*,*) "input graph: ", trim(inGraph)
      write(*,*) "output graph: ", trim(outGraph)
    end if
  end subroutine
end module

program main
  use engpar
  use splitHelperFns
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: ierr, self
  integer(ENGPAR_EDGE_T) :: edgeType
  type(c_ptr) :: graph, splitInput, diffusiveInput
  logical(C_BOOL) :: inSmall
  real(C_DOUBLE) :: tol, stepfactor
  integer :: smallSize, verbosity, newComm
  integer(ENGPAR_PART_T), dimension(1) :: ranks ! only used with LOCAL_PARMETIS
  character(len=256) :: inGraphFileName, outGraphFileName, splitMethod

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call cengpar_initialize()

  call getArgs(self, inGraphFileName, outGraphFileName, smallSize)

  inSmall = .false.
  call switchToOriginals(smallSize,inSmall,newComm)

  ! Switch the internal communicator (this changes PCU so use PCU_Comm_... with caution)
  call cengpar_setftncommunicator(newComm)

  graph = cengpar_createEmptyGraph()
 
  if (inSmall) then
    ! Only the original parts will construct the graph
    call cengpar_loadFromFile(graph,cleanString(inGraphFileName))
    call cengpar_evaluatePartition(graph);
  end if

  ! Create the split input
  edgeType = 0
  tol = 1.05
  ranks = 0
  splitInput = cengpar_createGlobalSplitInput(graph,newComm,MPI_COMM_WORLD,inSmall,tol,edgeType)

  ! Perform split 
  splitMethod = c_char_"GLOBAL_PARMETIS"//c_null_char
  call cengpar_split(splitInput,splitMethod);
  call cengpar_checkValidity(graph);
  if (self==0) write(*,*) 'After Split'
  call cengpar_evaluatePartition(graph)
  call cengpar_saveToFile(graph,cleanString(outGraphFileName))
  call MPI_Comm_free(newComm, ierr)

  call cengpar_destroyGraph(graph);
  call cengpar_finalize()
  call mpi_finalize(ierr)
  stop
end

