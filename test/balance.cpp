#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <cstring>
#include "buildGraphs.h"

bool cmpebin(char* str) {
  return strlen(str)>4&&strcmp(str+strlen(str)-5,".ebin")==0;             
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if ( argc!=1 &&argc!= 2 && argc!=3) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph>\n", argv[0]);
      printf("Usage: %s <bgd_prefix> [second edge type]\n", argv[0]);
      printf("Usage: %s\n",argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }
  agi::Ngraph* g=NULL;
  if (argc==1) {
    g= buildEmptyGraph();
  }
  else if (cmpebin(argv[1]))
    g= agi::createBinGraph(argv[1]);
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }

  if (g->hasCoords()) {
    std::string filename = "before";
    agi::writeVTK(g,filename.c_str());
  }
  engpar::evaluatePartition(g);

  double step_factor = 0.1;
  engpar::DiffusiveInput* input = engpar::createDiffusiveInput(g,step_factor);
  input->addPriority(0,1.1);
  if (argc==3) {
    input->addPriority(1,1.1);
  }
  input->addPriority(-1,1.1);

  input->maxIterationsPerType=50;
  input->maxIterations=75;
  //Create the balancer
  engpar::balance(input,2);

  //Evaluate the new partition
  engpar::evaluatePartition(g);

  //Ensure the graph is still valid
  agi::checkValidity(g);

  //Ensure the graph was balanced to the target tolerance
  assert(engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,-1))<1.2);
  assert(engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,0))<1.2);
  
  if (g->hasCoords()) {
    std::string filename = "after";
    agi::writeVTK(g,filename.c_str());
  }

  //Migration of original data structure
  if (argc>2) {
    agi::PartitionMap* map = g->getPartition();
    delete map;
  }
  
  //Destroy graph
  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
