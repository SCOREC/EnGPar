#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <PCU.h>
#include <sys/types.h>
#include <unistd.h>
void switchToOriginals(int smallSize, bool& isOriginal, MPI_Comm& newComm);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if (argc != 4) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph_file> <original_num_parts> <out_file>\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }

  //Application code:
  bool isOriginal = false;
  MPI_Comm newComm;
  int smallSize = atoi(argv[2]);
  switchToOriginals(smallSize, isOriginal,newComm);

  //Switch the internal communicator (this changes PCU so use PCU_Comm_... with caution)
  EnGPar_Switch_Comm(newComm);

  //Calls to EnGPar:
  agi::Ngraph* g = agi::createEmptyGraph();
  
  if (isOriginal) {
    //Only the original parts will construct the graph
    g->loadFromFile(argv[1]);
  }

  //Create the input
  double tolerance = 1.05;
  agi::etype t = 0;
  engpar::Input* input = engpar::createGlobalSplitInput(g,newComm,MPI_COMM_WORLD, isOriginal,
                                                        tolerance,t);

  if (isOriginal)
    engpar::evaluatePartition(g);
  engpar::split(input,engpar::GLOBAL_PARMETIS);
  engpar::evaluatePartition(g);

  //Application continues:
  MPI_Comm_free(&newComm);

  g->saveToFile(argv[3]);

  agi::destroyGraph(g);
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}


void switchToOriginals(int smallSize, bool& isOriginal, MPI_Comm& newComm) {
  int self = PCU_Comm_Self();
  int group;
  int groupRank;
  isOriginal = self<smallSize;

  if (isOriginal) {
    group=0;
    groupRank=self;
  }
  else {
    group = 1;
    groupRank = 0;
  }
  MPI_Comm_split(MPI_COMM_WORLD,group,groupRank,&newComm);
}
