#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <cstring>

bool cmpebin(char* str) {
  return strlen(str)>4&&strcmp(str+strlen(str)-5,".ebin")==0;             
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if ( argc!= 2 && argc!=3) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph>\n", argv[0]);
      printf("Usage: %s <bgd_prefix> [second edge type]\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }
  agi::Ngraph* g=NULL;
  if (cmpebin(argv[1]))
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

  engpar::Input* input = new engpar::Input(g);
  input->priorities.push_back(0);
  input->tolerances.push_back(1.1);
  if (argc==3) {
    input->priorities.push_back(1);
    input->tolerances.push_back(1.1);
  }
  input->priorities.push_back(-1);
  input->tolerances.push_back(1.1);

  input->step_factor=.1;

  input->useDistanceQueue=true;
  //Create the balancer
  agi::Balancer* balancer = engpar::makeBalancer(input);
  balancer->balance(1.1);

  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;

  //Ensure the graph is still valid
  agi::checkValidity(g);

  if (g->hasCoords()) {
    std::string filename = "after";
    agi::writeVTK(g,filename.c_str());
  }

  //Migration of original data structure
  if (argc>2) {
    g->getPartition();

  }
  
  //Destroy graph
  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
