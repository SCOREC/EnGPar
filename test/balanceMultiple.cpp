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
  
  if (argc!= 2) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <bgd_prefix>\n", argv[0]);
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

  engpar::evaluatePartition(g);
  double step_factor = 0.1;
  engpar::Input* input = engpar::createDiffusiveInput(g,step_factor);
  input->addPriority(-1,1.1);

  //Create the balancer
  agi::Balancer* balancer = engpar::makeBalancer(input,1);
  balancer->balance(1.1);

  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;

  //Ensure the graph is still valid
  agi::checkValidity(g);

  //Alter weights so odd processes are heavier than even processes
  agi::wgt_t* weights = new agi::wgt_t[g->numLocalVtxs()];
  agi::wgt_t w = PCU_Comm_Self()%2+1;

  for (agi::lid_t i =0;i<g->numLocalVtxs();i++) {
    weights[i] = w;
  }
  g->setWeights(weights);
  delete [] weights;

  g->resetOwnership();
  
  engpar::evaluatePartition(g);
  input = engpar::createDiffusiveInput(g,step_factor);
  input->addPriority(-1,1.1);


  //Create the balancer
  balancer = engpar::makeBalancer(input,1);
  balancer->balance(1.1);

  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;

  //Ensure the graph is still valid
  agi::checkValidity(g);

  
  //Destroy graph
  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
