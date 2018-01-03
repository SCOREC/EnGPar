#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <engpar_split.h>
#include <PCU.h>
#include <sys/types.h>
#include <unistd.h>
void switchToOriginals(int split_factor, bool& isOriginal, MPI_Comm& newComm);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if (argc != 4) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph_file> <split_factor> <out_file>\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }

  //Application code:
  bool isOriginal = false;
  MPI_Comm newComm;
  int split_factor = atoi(argv[2]);
  switchToOriginals(split_factor, isOriginal,newComm);

  //Calls to EnGPar:
  agi::Ngraph* g = agi::createEmptyGraph();

  //Create the input (this sets up the internal communicators,
  //                    so this must be done before graph construction.)
  double tolerance = 1.05;
  agi::etype t = 0;
  engpar::Input* input_s = engpar::createSplitInput(g,newComm,MPI_COMM_WORLD, isOriginal,
                                                  split_factor,tolerance,t);
  
  if (isOriginal) {
    //Only the original parts will construct the graph
    g->loadFromFile(argv[1]);
    if (!PCU_Comm_Self())
      printf("\nBefore Split\n");
    engpar::evaluatePartition(g);
  }

  engpar::split(input_s,engpar::GLOBAL_PARMETIS);
  if (!PCU_Comm_Self())
    printf("\nAfter Split\n");
  engpar::evaluatePartition(g);
  delete input_s;
  
  //Create the input for diffusive load balancing (vtx>element)
  double step_factor = 0.1;
  engpar::Input* input_d = engpar::createDiffusiveInput(g,step_factor);
  input_d->addPriority(0,1.03);
  input_d->addPriority(-1,1.03);
  if (!PCU_Comm_Self())
    printf("\n");

  //Create and run the balancer
  agi::Balancer* balancer = engpar::makeBalancer(input_d,0);
  balancer->balance(1.1);

  if (!PCU_Comm_Self())
    printf("\nAfter Balancing\n");
  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;

  //Ensure the graph is still valid
  agi::checkValidity(g);
  
  //Application continues:
  MPI_Comm_free(&newComm);

  g->saveToFile(argv[3]);

  agi::destroyGraph(g);
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}


void switchToOriginals(int split_factor, bool& isOriginal, MPI_Comm& newComm) {
  int self = PCU_Comm_Self();
  int group;
  int groupRank;
  isOriginal = self%split_factor==0;

  if (isOriginal) {
    group=0;
    groupRank=self/split_factor;
  }
  else {
    group = 1;
    groupRank = 0;
  }
  MPI_Comm_split(MPI_COMM_WORLD,group,groupRank,&newComm);
}
