#include <binGraph.h>
#include <mpi.h>
#include <PCU.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <ZoltanCutVertex.h>

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2&&argc!=3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file> [vertex_partition_file]",argv[0]);
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }
  agi::binGraph* g;
  if (argc==2)
    g = new agi::binGraph(argv[1]);
  else
    g = new agi::binGraph(argv[1],argv[2]);
  zagi::ZoltanCutVertex* ptn = new zagi::ZoltanCutVertex(g,2);
  ptn->run();

  delete ptn;
  delete g;

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  PCU_Comm_Free();
  MPI_Finalize();

}
