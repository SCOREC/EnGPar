#include <ngraph.h>
#include <binGraph.h>
#include <PCU.h>
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file>",argv[0]);
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }

  //Construct Graph from binary file
  agi::Ngraph* g;

  g = agi::createBinGraph(argv[1]);

  //Run BFS algorithsm here
  
  destroyGraph(g);

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}
