#include <Kokkos_Core.hpp>
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
  KOKKOS_FOR_VERTS(g,v) {
    printf("%lu\n",g->globalID(v));
    printf("%lu\n",g->localID(v));
  }
  KOKKOS_END_FOR()
  
  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
