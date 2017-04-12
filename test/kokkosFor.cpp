#include <Kokkos_Core.hpp>
#include <binGraph.h>
#include <mpi.h>
#include <PCU.h>
#include <stdio.h>


int main(int argc, char* argv[]) {

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file>",argv[0]);
    PCU_Comm_Free();
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
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;

}
