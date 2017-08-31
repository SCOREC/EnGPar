#include <Kokkos_Core.hpp>
#include <binGraph.h>
#include <mpi.h>
#include <stdio.h>
#include <engpar_support.h>
#include <ngraph.h>
#include <pngraph.h>

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

  typedef Kokkos::View< agi::wgt_t* > wgt_view;
  wgt_view x_dev("x", png->num_local_verts);

  Kokkos::parallel_for(png->num_local_verts,
      KOKKOS_LAMBDA(int i) {
        x_dev(i) = -1;
      }
  );

  // create a device array with the local weights
  wgt_view weights_dev("local_weights", png->num_local_verts);
  wgt_view::HostMirror weights_host = Kokkos::create_mirror_view( weights_dev );
  // initialize it on the host
  for(int i=0; i<png->num_local_verts; i++) {
    weights_host(i) = png->local_weights[i];
  }
  // deep copy from the host to the device
  Kokkos::deep_copy( weights_dev, weights_host );

  // write the weights into the x array on the device
  Kokkos::parallel_for(png->num_local_verts,
      KOKKOS_LAMBDA(int i) {
        x_dev(i) = weights_dev(i);
      }
  );

  // deep copy the x array from the device to the host
  wgt_view::HostMirror x_host = Kokkos::create_mirror_view( x_dev );
  Kokkos::deep_copy( x_host, x_dev );

  // check that the copy was successful
  for(int i=0; i<png->num_local_verts; i++) {
    assert(x_host[i] == png->local_weights[i]);
  }

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
