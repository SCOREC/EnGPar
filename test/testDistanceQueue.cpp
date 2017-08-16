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
    if (PCU_Comm_Peers()==2)
      g = buildDisconnected2Graph();
    else
      g=buildHyperGraphLine();
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }
  PCU_Barrier();

  engpar::Queue* q = engpar::createDistanceQueue(g);
  for (unsigned int i=0;i<q->size();i++)
    printf("%d %ld\n",PCU_Comm_Self(),g->globalID((*q)[i]));
  delete q;

  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n"); 

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}
