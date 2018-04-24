!---------------------------------------------------------------------------
!> @file engpar.f90
!> @brief EnGPar FORTRAN interface using iso_c_binding
!---------------------------------------------------------------------------
module engpar
  use :: iso_c_binding
#include "agi_types.h"
  public
  integer, parameter :: ENGPAR_GID_T  = AGI_GID_FT
  integer, parameter :: ENGPAR_LID_T  = AGI_LID_FT
  integer, parameter :: ENGPAR_WGT_T  = AGI_WGT_FT
  integer, parameter :: ENGPAR_PART_T = AGI_PART_FT
  integer, parameter :: ENGPAR_EDGE_T = AGI_EDGE_FT
  interface
  !---------------------------------------------------------------------------
  !> @brief initialize engar, call this before any other engpar api
  !---------------------------------------------------------------------------
  subroutine cengpar_initialize() bind(C, NAME='cengpar_initialize')
    use :: iso_c_binding
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief finalize engar, no engpar apis may be called after this
  !---------------------------------------------------------------------------
  subroutine cengpar_finalize() bind(C, NAME='cengpar_finalize')
    use :: iso_c_binding
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the MPI communicator used by engpar
  !> @param comm(in) MPI communicator
  !---------------------------------------------------------------------------
  subroutine cengpar_setftncommunicator(comm) &
             bind(C, NAME='cengpar_setftncommunicator')
    use :: iso_c_binding
    integer(c_int), value :: comm
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief create an empty engpar graph
  !> @return empty graph
  !---------------------------------------------------------------------------
  function cengpar_createEmptyGraph() &
             bind(C, NAME='cengpar_createEmptyGraph')
    use :: iso_c_binding
    type(c_ptr) cengpar_createEmptyGraph
  end function

  !---------------------------------------------------------------------------
  !> @brief create input for splitting
  !> @param graph(in) engpar graph
  !> @param smallComm(in) small mpi communicator
  !> @param largeComm(in) large mpi communicator
  !> @param isOrig(in) true if process is part of the small communicator
  !> @param splitFactor(in) = |large comm| / |small comm|
  !> @param tol(in) target imbalance tolerance used for splitting
  !> @param ranks(in) list of MPI ranks to use for local splitting
  !> @param edgeType(in) mesh entity dimension used to create graph edges
  !> @return input for splitting
  !---------------------------------------------------------------------------
  function cengpar_createLocalSplitInput(graph,smallComm,largeComm,isOrig,splitFactor, &
             tol,ranks,edgeType) &
             bind(C, NAME='cengpar_createLocalSplitInput')
    use :: iso_c_binding
    type(c_ptr) :: cengpar_createLocalSplitInput
    type(c_ptr), value :: graph
    integer(c_int), value :: smallComm
    integer(c_int), value :: largeComm
    logical(c_bool), intent(in), value :: isOrig
    integer(c_int), value :: splitFactor
    real(c_double), value :: tol
    integer(AGI_PART_FT), intent(in), dimension(splitFactor) :: ranks
    integer(AGI_EDGE_FT), value :: edgeType
  end function
  !---------------------------------------------------------------------------
  !> @brief create input for splitting
  !> @param graph(in) engpar graph
  !> @param smallComm(in) small mpi communicator
  !> @param largeComm(in) large mpi communicator
  !> @param isOrig(in) true if process is part of the small communicator
  !> @param tol(in) target imbalance tolerance used for splitting
  !> @param edgeType(in) mesh entity dimension used to create graph edges
  !> @param ranks(in) list of MPI ranks to use for local splitting
  !> @return input for splitting
  !---------------------------------------------------------------------------
  function cengpar_createGlobalSplitInput(graph,smallComm,largeComm,isOrig, &
             tol,edgeType) &
             bind(C, NAME='cengpar_createGlobalSplitInput')
    use :: iso_c_binding
    type(c_ptr) :: cengpar_createGlobalSplitInput
    type(c_ptr), value :: graph
    integer(c_int), value :: smallComm
    integer(c_int), value :: largeComm
    logical(c_bool), intent(in), value :: isOrig
    real(c_double), value :: tol
    integer(AGI_EDGE_FT), value :: edgeType
  end function

  !---------------------------------------------------------------------------
  !> @brief create input for the diffusive balaancer
  !> @param graph(in) engpar graph
  !> @param stepfactor(in) controls how much weight is migrated in each
  !>                       iteration, start with a setting of 0.1
  !> @return input for balancing
  !---------------------------------------------------------------------------
  function cengpar_createDiffusiveInput(graph,stepfactor) &
             bind(C, NAME='cengpar_createDiffusiveInput')
    use :: iso_c_binding
    type(c_ptr) :: cengpar_createDiffusiveInput
    type(c_ptr), value :: graph
    real(c_double), value :: stepfactor
  end function
  !---------------------------------------------------------------------------
  !>  @brief add entity type to balance
  !>  @remark add edge type t with tolerance tol to the priority list, 
  !>  if edge type t appears edge type t' then t has a higher priority than t'
  !>  @param input (in/out) diffusive input created with 'cengpar_createDiffusiveInput'
  !>  @param etype (in) -1 represents graph vertices, 0-MAX_TYPES represent edge types
  !>  @param tol (in) 1.0 + alpha allows for an imbalance of alpha percent,
  !>                  alpha = 0 requests perfect balancing
  !---------------------------------------------------------------------------
  subroutine cengpar_addPriority(input,etype,tol) &
             bind(C, NAME='cengpar_addPriority')
    use :: iso_c_binding
    type(c_ptr), value :: input
    integer(c_int), value :: etype
    real(c_double), value :: tol
  end subroutine
  !---------------------------------------------------------------------------
  !>  @brief balance the graph based on the specified input
  !>  @remark the engparInput object is destroyed
  !>  @param input (in) the object containing diffusive balancer inputs
  !>  @param verbosity (In) as verbosity increases the amount of output
  !>                        increases, min verbosity=0, max verbosity=3
  !---------------------------------------------------------------------------
  subroutine cengpar_balance(input,verbosity) &
             bind(C, NAME='cengpar_balance')
    use :: iso_c_binding
    type(c_ptr), value :: input
    integer(C_INT), intent(in), value :: verbosity
  end subroutine

  !---------------------------------------------------------------------------
  !> @brief load graph from file
  !> @param graph(in/out) empty engpar graph
  !> @param fileName (in) path to graph file to load (prefix only, no '_#.bgd')
  !---------------------------------------------------------------------------
  subroutine cengpar_loadFromFile(graph,fileName) &
             bind(C, NAME='cengpar_loadFromFile')
    use :: iso_c_binding
    type(c_ptr), value :: graph
    character(c_char), intent(in) :: fileName(*)
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief save graph to file
  !> @param graph(in) engpar graph
  !> @param fileName (in) path to save graph file (prefix only, no '_#.bgd')
  !---------------------------------------------------------------------------
  subroutine cengpar_saveToFile(graph,fileName) &
             bind(C, NAME='cengpar_saveToFile')
    use :: iso_c_binding
    type(c_ptr), value :: graph
    character(c_char), intent(in) :: fileName(*)
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief print partition info
  !> @param graph(in) engpar graph
  !---------------------------------------------------------------------------
  subroutine cengpar_evaluatePartition(graph) &
             bind(C, NAME='cengpar_evaluatePartition')
    use :: iso_c_binding
    type(c_ptr), value :: graph
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief split the graph
  !> @param input (in/out) split input created with 'cengpar_create*SplitInput'
  !> @param splitMethod (in) string matching one of the engpar_split.h enums
  !---------------------------------------------------------------------------
  subroutine cengpar_split(input,splitMethod) &
             bind(C, NAME='cengpar_split')
    use :: iso_c_binding
    type(c_ptr), value :: input
    character(c_char), intent(in) :: splitMethod(*)
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief construct vertices
  !> @param graph(in/out) empty engpar graph
  !> @param isHg(in) true if hypergraph, false otherwise
  !> @param verts(in) array of graph vertex ids
  !> @param weights(in) array of graph vertex weights
  !> @param nverts(in) number of graph vertices
  !---------------------------------------------------------------------------
  subroutine cengpar_constructVerts(graph,isHg,verts,weights,nverts) &
             bind(C, NAME='cengpar_constructVerts')
    use :: iso_c_binding
    type(c_ptr), value :: graph
    logical(c_bool), intent(in), value :: isHg
    integer(AGI_GID_FT), intent(in), dimension(nverts) :: verts
    real(AGI_WGT_FT), intent(in), dimension(nverts) :: weights
    integer(C_INT), intent(in), value :: nverts
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief construct edges
  !> @param graph(in/out) engpar graph
  !> @param edges(in) array of graph edge ids
  !> @param degs(in) array of graph edge degrees (number of vertices bounding each edge)
  !> @param weights(in) array of graph edge weights
  !> @param pins(in) array of graph vertices bounding each edge
  !> @param nedges(in) number of graph edges
  !> @param npins(in) number of pins
  !---------------------------------------------------------------------------
  subroutine cengpar_constructEdges(graph,edges,degs,weights,pins,nedges,npins) &
             bind(C, NAME='cengpar_constructEdges')
    use :: iso_c_binding
    type(c_ptr), value :: graph
    integer(AGI_GID_FT), intent(in), dimension(nedges) :: edges
    integer(AGI_LID_FT), intent(in), dimension(nedges) :: degs
    real(AGI_WGT_FT), intent(in), dimension(nedges) :: weights
    integer(AGI_GID_FT), intent(in), dimension(npins) :: pins
    integer(C_INT), intent(in), value :: nedges
    integer(C_INT), intent(in), value :: npins
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief construct ghost vertices
  !> @remark for each $e(u,v)$, a ghost vertex is needed if the owner of $u$ is not the owner of $v$.
  !>  This must be called after edges are created.
  !> @param graph(in/out) engpar graph
  !> @param verts(in) array of ghosted graph vertices (global ids)
  !> @param owners(in) array of ghost graph vertex owners; the part id that owns each graph vertex
  !> @param nghosts(in) number of ghost vertices
  !---------------------------------------------------------------------------
  subroutine cengpar_constructGhosts(graph,verts,owners,nghosts) &
             bind(C, NAME='cengpar_constructGhosts')
    use :: iso_c_binding
    type(c_ptr), value :: graph
    integer(AGI_GID_FT), intent(in), dimension(nghosts) :: verts
    integer(AGI_PART_FT), intent(in), dimension(nghosts) :: owners
    integer(C_INT), intent(in), value :: nghosts
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief check the validity of the graph
  !> @param graph(in/out) engpar graph
  !---------------------------------------------------------------------------
  subroutine cengpar_checkValidity(graph) &
             bind(C, NAME='cengpar_checkValidity')
    use :: iso_c_binding
    type(c_ptr), value :: graph
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief destroy the graph
  !> @param graph(in/out)  engpar graph
  !---------------------------------------------------------------------------
  subroutine cengpar_destroyGraph(graph) &
             bind(C, NAME='cengpar_destroyGraph')
    use :: iso_c_binding
    type(c_ptr), value :: graph
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief balance the graph vertices
  !> @remark the imbalance is specified as the maximum vertex weight across all parts
  !>          divided by the average vertex weight across all parts
  !> @param graph(in/out)  engpar graph
  !> @param tol(in) target maximum imbalance max(vertex weight) / average(vertex weight)
  !> @param stepfactor(in) controls how much weight is migrated in each
  !>                       iteration, start with a setting of 0.1
  !> @param verbosity(in) the level of output; the higher the value the more output
  !---------------------------------------------------------------------------
  subroutine cengpar_balanceVertices(graph,tol,stepfactor,verbosity) &
             bind(C, NAME='cengpar_balanceVertices')
    use :: iso_c_binding
    type(c_ptr), value :: graph
    real(C_DOUBLE), intent(in), value :: tol, stepfactor
    integer(C_INT), intent(in), value :: verbosity
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the migration schedule computed by a balancing method
  !> @param graph(in) engpar graph
  !> @param verts(in/out) array of graph vertex global ids
  !> @param parts(in/out) array of destination part ids for each vertex
  !> @param nverts(in) number of graph vertices
  !---------------------------------------------------------------------------
  subroutine cengpar_getPartition(graph,verts,parts,nverts) &
             bind(C, NAME='cengpar_getPartition')
    use :: iso_c_binding
    type(c_ptr), value :: graph
    integer(AGI_GID_FT), intent(in), dimension(nverts) :: verts
    integer(AGI_PART_FT), intent(in), dimension(nverts) :: parts
    integer(C_INT), intent(in), value :: nverts
  end subroutine
  end interface
end module
