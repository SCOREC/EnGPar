#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <cstring>
#include "buildGraphs.h"
#include <engpar_diffusive_input.h>
#include <engpar_metrics.h>

bool cmpebin(char* str) {
  return strlen(str)>4&&strcmp(str+strlen(str)-5,".ebin")==0;             
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if (argc < 3 || argc > 5) {
    if (!PCU_Comm_Self()) {
      printf("Usage: %s <graph.ebin> <tolerance> [verbosity] [save_prefix]\n", argv[0]);
      printf("Usage: %s <bgd_prefix> <tolerance> [verbosity] [save_prefix]\n", argv[0]);
    }
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }
  agi::Ngraph* g=NULL;
  if (cmpebin(argv[1]))
    g= agi::createBinGraph(argv[1]);
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }

  double tol = atof(argv[2]);
  
  //Evaluate the new partition
  engpar::evaluatePartition(g);

  double step_factor = 0.1;
  engpar::DiffusiveInput* input = engpar::createDiffusiveInput(g,step_factor);
  for (agi::etype t = 0; t < g->numEdgeTypes(); ++t)
    input->addPriority(t,tol);
  input->addPriority(-1,tol);

  int verbosity = 1;
  if (argc > 3)
    verbosity = atoi(argv[3]);
  //Create the balancer
  engpar::balance(input,verbosity);

  if (argc > 4) {
    g->saveToFile(argv[4]);
  }
  //Evaluate the new partition
  engpar::evaluatePartition(g);
  
  //Destroy graph
  agi::destroyGraph(g);

  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}
