#include <engpar_support.h>
#include <binGraph.h>
#include <engpar.h>
#include <cstring>

bool cmpebin(char* str) {
  return strlen(str)>4&&strcmp(str+strlen(str)-5,".ebin")==0;             
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if ( argc!=3 &&argc != 4) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph> <step factor> [verbosity]\n", argv[0]);
      printf("Usage: %s <bgd_prefix> <step factor> [verbosity]\n", argv[0]);
    }
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  //Load mesh
  agi::Ngraph* g= NULL;
  double step_factor = atof(argv[2]);
  int verbosity=0;
  if (argc>3)
    verbosity = atoi(argv[3]);
  if (cmpebin(argv[1])) {
    g = agi::createBinGraph(argv[1]);
  }
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }
  
  engpar::evaluatePartition(g);
  
  //Create the balancer
  agi::Balancer* balancer = engpar::makeVtxBalancer(g,step_factor,verbosity);
  balancer->balance(1.01);
  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;
  
  //Destroy graph
  agi::destroyGraph(g);
  
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
