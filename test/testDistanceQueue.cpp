#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include "buildGraphs.h"
#include "../partition/src/engpar_queue.h"
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  agi::Ngraph* g;
  if (argc==1)
    g = buildHyperGraphLine();
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }
  PCU_Barrier();

  engpar::Queue* q = engpar::createDistanceQueue(g);
  for (unsigned int i=0;i<q->size();i++)
    printf("%lu\n",g->localID((*q)[i]));
  delete q;

  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n"); 

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}
