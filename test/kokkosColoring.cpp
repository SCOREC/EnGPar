#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <binGraph.h>
#include <mpi.h>
#include <stdio.h>
#include <engpar_support.h>

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

  agi::PNgraph* pg = g->publicize();

  agi::lid_t numverts = pg->num_local_verts;

  agi::lid_t numedges = pg->num_local_edges[0];  // Add support for multiple edge types

  agi::lid_t* degree_list = pg->degree_list[0]; // ^

  agi::lid_t* edge_list = pg->edge_list[0]; // ^

  // Create views for each part of the graph data-structure
  
  Kokkos::View<agi::lid_t*> degree_view ("degree_view",numverts);

  Kokkos::View<agi::lid_t*> edge_view ("edge_view",numedges);
 
  // Use parllel loops to fill the views
  Kokkos::parallel_for(numverts+1, KOKKOS_LAMBDA( const int i ) {
    degree_view(i) = degree_list[i];
  });

  Kokkos::parallel_for(numedges, KOKKOS_LAMBDA( const int i) {
    edge_view(i) = edge_list[i];
  });

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
