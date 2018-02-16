#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include "buildGraphs.h"
#include "sellCSigma.h"
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

  agi::lid_t C = 4;
  agi::lid_t sigma = 64;
  agi::Ngraph* scg = ssg::convertFromAGI(g,C,sigma);

  agi::destroyGraph(g);

  agi::destroyGraph(scg);

  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n"); 

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}
