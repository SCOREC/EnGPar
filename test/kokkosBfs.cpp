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
  fprintf(stderr, "num_types %d\n", png->num_types);
  const int edgeType = 0;
  assert(png->isHyperGraph);

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

  // copy the vertex degree list to the device
  int_array vtx_to_net_offsets("vtx to net offsets", png->num_local_verts+1);
  int_array::HostMirror vtx_to_net_offsets_h = Kokkos::create_mirror_view( vtx_to_net_offsets );
  for(int i=0; i<png->num_local_verts+1; i++)
    vtx_to_net_offsets_h(i) = png->degree_list[edgeType][i];
  Kokkos::deep_copy( vtx_to_net_offsets, vtx_to_net_offsets_h );

  // copy the adjacent vertex list to the device
  int_array vtx_to_net("vtx to net adjacencies", png->num_local_edges[edgeType]);
  int_array::HostMirror vtx_to_net_h = Kokkos::create_mirror_view( vtx_to_net );
  for(int i=0; i<png->num_local_edges[edgeType]; i++)
    vtx_to_net_h(i) = png->edge_list[edgeType][i];
  Kokkos::deep_copy( vtx_to_net, vtx_to_net_h );

  // create the queues (device work arrays)
  int_array queue("queue", png->num_local_verts*QUEUE_MULTIPLIER);
  int_array queue_next("queue next", png->num_local_verts*QUEUE_MULTIPLIER);

  // run the traversal
  fwbw_baseline(vtx_to_net_offsets, vtx_to_net, root, visited, queue, queue_next);

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
