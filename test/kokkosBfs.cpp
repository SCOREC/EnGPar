#include <Kokkos_Core.hpp>
#include <binGraph.h>
#include <mpi.h>
#include <stdio.h>
#include <engpar_support.h>
#include <ngraph.h>
#include <pngraph.h>
#include "kokkosBfs.h"
#include "bfs_baseline.cpp"

bool cmpebin(char* str) {
  return strlen(str)>4&&strcmp(str+strlen(str)-5,".ebin")==0;
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

  agi::Ngraph* g=NULL;
  if (cmpebin(argv[1]))
    g = agi::createBinGraph(argv[1]);
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }

  agi::PNgraph* png = g->publicize();
  fprintf(stderr, "num_local_verts %d\n", png->num_local_verts);

  typedef Kokkos::View< agi::lid_t > lid_type;

  int_type num_local_verts("num_local_verts", 1);
  int_type::HostMirror num_local_verts_h = Kokkos::create_mirror_view( num_local_verts );
  // set the host value
  num_local_verts_h() = png->num_local_verts;
  // copy to device
  Kokkos::deep_copy( num_local_verts, num_local_verts_h );

  bool_array visited("visited", png->num_local_verts);

  int_type root("root", 0); //pick the first vertex for now

  // initialize visited
  fwbw_init(num_local_verts, visited, root);

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
