#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include "buildGraphs.h"
#include "../partition/src/engpar_queue.h"
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Debug_Open();
  
  agi::Ngraph* g = buildHyperGraphLine();
  PCU_Barrier();

  engpar::Queue* q = engpar::createDistanceQueue(g);

  delete q;

  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n"); 

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}
