#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include "buildGraphs.h"
#include "../partition/Diffusive/engpar_diffusive_input.h"
#include "../partition/Diffusive/src/engpar_queue.h"
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

  engpar::Input* input = engpar::createDiffusiveInput(g,0);
  engpar::DiffusiveInput* inp = static_cast<engpar::DiffusiveInput*>(input);
  engpar::Queue* q = engpar::createDistanceQueue(inp);
  if (!PCU_Comm_Self()) {
    for (unsigned int i=0;i<q->size();i++)
      printf("%d %ld\n",PCU_Comm_Self(),g->globalID((*q)[i]));
  }
  delete q;
  delete input;
  
  if (g->hasCoords()) {
    std::string filename = "graph";
    agi::writeVTK(g,filename.c_str());
  }

  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n"); 

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}
