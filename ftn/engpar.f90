!---------------------------------------------------------------------------
!> @file engpar.f90
!> @brief EnGPar FORTRAN interface using iso_c_binding
!---------------------------------------------------------------------------
module engpar
  use :: iso_c_binding
#include "agi_types.h"
  public
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
  !> @param stepFactor(in) controls how much weight is migrated in each
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
  !> @param verts(in) array of graph vertex global ids
  !> @param parts(in) array of destination part ids for each vertex
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
