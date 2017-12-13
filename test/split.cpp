#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <engpar_split.h>
#include <PCU.h>
#include <sys/types.h>
#include <unistd.h>
bool switchToOriginals(int split_factor);
void switchToAll(agi::Ngraph* g, bool isOriginal);

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

  bool isOriginal = switchToOriginals(atoi(argv[2]));
  agi::Ngraph* g = agi::createEmptyGraph();
  agi::Migration* plan;
  if (isOriginal) {
    g->loadFromFile(argv[1]);
    engpar::evaluatePartition(g);
    plan =engpar::split(g,atoi(argv[2]),engpar::GLOBAL_PARMETIS);
  }
  else
    plan = new agi::Migration(g);
  switchToAll(g,isOriginal);
  g->migrate(plan);
  engpar::evaluatePartition(g);

  g->saveToFile(argv[3]);
  
  agi::destroyGraph(g);
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}


bool switchToOriginals(int split_factor) {
  int self = PCU_Comm_Self();
  int group;
  int groupRank;
  bool isOriginal = self%split_factor==0;

  if (isOriginal) {
    group=0;
    groupRank=self/split_factor;
  }
  else {
    group = 1;
    groupRank = 0;
  }
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD,group,groupRank,&groupComm);
  PCU_Switch_Comm(groupComm);
  return isOriginal;
}

void switchToAll(agi::Ngraph* g, bool isOriginal)
{
  MPI_Comm prevComm = PCU_Get_Comm();

  //Expands the view of the graph from s parts to t parts
  if (isOriginal)
    engpar::expandParts(g,MPI_COMM_WORLD);
  else
    PCU_Switch_Comm(MPI_COMM_WORLD);

  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}
