#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosGraph_graph_color.hpp>
#include <KokkosKernels_Handle.hpp>
#include <binGraph.h>
#include <mpi.h>
#include <stdio.h>
#include <engpar_support.h>
#include <cstring>

#include <iostream>

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

  agi::Ngraph* g = agi::createEmptyGraph();
  g->loadFromFile(argv[1]);

  agi::PNgraph* pg = g->publicize();

  const agi::lid_t numverts = pg->num_local_verts;

  const agi::lid_t numedges = pg->num_local_edges[0];  // Add support for multiple edge types

  const agi::lid_t* degree_list = pg->degree_list[0]; // ^

  const agi::lid_t* edge_list = pg->edge_list[0]; // ^

  // Create views for each part of the graph data-structure
  
  Kokkos::View<agi::lid_t*> degree_view ("degree_view",numverts+1);

  Kokkos::View<agi::lid_t*> edge_view ("edge_view",numedges);
 
  Kokkos::View<agi::lid_t*> vals_view("vals_view",numedges);

  // Use parllel loops to fill the views
  
  Kokkos::parallel_for(numverts+1, KOKKOS_LAMBDA( const int i ) {
    degree_view(i) = degree_list[i];
  });

  Kokkos::parallel_for(numedges, KOKKOS_LAMBDA( const int i) {
    edge_view(i) = edge_list[i];
    vals_view(i) = 1;
  });

  KokkosSparse::CrsMatrix<agi::lid_t, agi::lid_t, Kokkos::Serial::device_type, void, int>("graph", numverts, numverts, numedges, vals_view, degree_view, edge_view);

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<agi::lid_t, agi::lid_t, agi::lid_t,
                                                           Kokkos::Serial::execution_space, 
                                                           Kokkos::Serial::memory_space, 
                                                           Kokkos::Serial::memory_space> KernelHandle;


  KernelHandle *kh = new KernelHandle();
  kh->set_team_work_size(16);
  kh->set_dynamic_scheduling(true);
  kh->create_graph_coloring_handle(KokkosGraph::COLORING_SERIAL);

  KokkosGraph::Experimental::graph_color_symbolic
	<KernelHandle, Kokkos::View<agi::lid_t*>, Kokkos::View<agi::lid_t*> >
	(kh, numverts, numverts, degree_view, edge_view);

  Kokkos::View<int*> vert_colors = kh->get_graph_coloring_handle()->get_vertex_colors();

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
  while ((e=g->iterate(eitr))) {
    int u = g->getIntTag(tag, g->u(e));
    int v = g->getIntTag(tag, g->v(e));
    assert(u!=v);
  }

  // Write the vtk files
//  std::string filename = "kokkos_color";
//  agi::writeVTK(g,filename.c_str(),tag,-1);

  // Finalize & Delete
  kh->destroy_graph_coloring_handle();
  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
