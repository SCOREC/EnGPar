#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
<<<<<<< HEAD
=======
#include <KokkosGraph_graph_color.hpp>
#include <KokkosKernels_Handle.hpp>
>>>>>>> lucas_coloring
#include <binGraph.h>
#include <mpi.h>
#include <stdio.h>
#include <engpar_support.h>
#include <cstring>

#include <iostream>


bool Check_Directed(agi::Ngraph* g, agi::etype t=0) {

  agi::PNgraph* pg = g->publicize();

  const agi::lid_t numverts = pg->num_local_verts; 

  const agi::lid_t* degree_list = pg->degree_list[t];

  const agi::lid_t* edge_list = pg->edge_list[t];

  for (int i=0; i<numverts; ++i) { 
    for (int j=degree_list[i]; j<degree_list[i+1]; ++j) {
      bool back_edge = false; 
      for (int k=degree_list[edge_list[j]]; k<degree_list[1+edge_list[j]]; ++k) {
        if (i == edge_list[k])
          back_edge = true;
      }
      if (!back_edge) {
        std::cout << "Directed graphs are not supported" << std::endl;
        return false;
      }
    }
  }
  return true;
}


void Color_Graph(agi::Ngraph* g, agi::etype t=0) {

  agi::PNgraph* pg = g->publicize();

  const agi::lid_t numverts = pg->num_local_verts;

  const agi::lid_t numedges = pg->num_local_edges[t];

  const agi::lid_t* degree_list = pg->degree_list[t];

  const agi::lid_t* edge_list = pg->edge_list[t];

  // Create views for each part of the graph data-structure
  
  Kokkos::View<agi::lid_t*> degree_view ("degree_view",numverts+1);

  Kokkos::View<agi::lid_t*> edge_view ("edge_view",numedges);

  // Use parllel loops to fill the views
  
  Kokkos::parallel_for(numverts+1, KOKKOS_LAMBDA( const int i ) {
    degree_view(i) = degree_list[i];
  });

  Kokkos::parallel_for(numedges, KOKKOS_LAMBDA( const int i) {
    edge_view(i) = edge_list[i];
  });

  typedef Kokkos::DefaultExecutionSpace exe_space;
  typedef KokkosSparse::CrsMatrix<agi::lid_t, agi::lid_t, Kokkos::Device<exe_space,exe_space::memory_space>, void, int> crsMat_t;
  typedef crsMat_t::StaticCrsGraphType graph_t;
  typedef graph_t::entries_type::non_const_type color_view_t;
  typedef graph_t::row_map_type lno_view_t;
  typedef graph_t::entries_type lno_nnz_view_t;
  typedef graph_t::entries_type::non_const_type  color_view_t;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle <agi::lid_t, agi::lid_t, agi::lid_t, 
            exe_space, exe_space::memory_space, exe_space::memory_space> KernelHandle; 

  // Create kernel handle
  KernelHandle* kh = new KernelHandle();
  kh->set_team_work_size(16);
  kh->set_dynamic_scheduling(true);
  kh->create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);

  // Run kokkos coloring
  KokkosGraph::Experimental::graph_color
	<KernelHandle, Kokkos::View<agi::lid_t*>, Kokkos::View<agi::lid_t*> >
	(kh, numverts, numverts, degree_view, edge_view);

  color_view_t vert_colors = kh->get_graph_coloring_handle()->get_vertex_colors();

  // ########## IMPL CHECK ##########
  assert( (0==KokkosKernels::Impl::kk_is_d1_coloring_valid <lno_view_t, lno_nnz_view_t, color_view_t, exe_space>
                (numverts, numverts, degree_view, edge_view, vert_colors)) );

  // Assign colors from kokkos graph to EnGPar graph
  agi::checkValidity(g);
  agi::GraphTag* tag = g->createIntTag(-1);
  agi::GraphVertex* v;
  agi::VertexIterator* vitr = g->begin();

  // Assigning colors to the original graph
  int i=0;
  while ((v=g->iterate(vitr))) {
    g->setIntTag(tag,v,vert_colors(i++));
  }
  
  // Check that coloring is valid
  agi::GraphEdge* e;
  agi::EdgeIterator* eitr = g->begin(0);
  int conflicts = 0;
  while ((e=g->iterate(eitr))) {
    int u = g->getIntTag(tag, g->u(e));
    int v = g->getIntTag(tag, g->v(e));
    if (u==v) ++conflicts;
  }
  g->destroy(eitr);
  assert(conflicts == 0);

  kh->destroy_graph_coloring_handle();
  
}

int main(int argc, char* argv[]) {

  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  if ( argc != 2 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file>",argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  Kokkos::initialize(argc,argv);

  agi::Ngraph* g = agi::createBinGraph(argv[1]);

  assert(Check_Directed(g));

  Color_Graph(g);

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
