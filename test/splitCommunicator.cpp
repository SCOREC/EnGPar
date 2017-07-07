#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include "buildGraphs.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  
  MPI_Comm new_comm;
  int new_size = PCU_Comm_Peers()/2;
  int old_self = PCU_Comm_Self();
  MPI_Comm_split(MPI_COMM_WORLD,old_self/new_size,
                 old_self%new_size,&new_comm);

  PCU_Switch_Comm(new_comm);

  printf("Self: %d\n",PCU_Comm_Self());

  if (old_self<new_size) {
    agi::Ngraph* g = buildGraph();

    agi::VertexIterator* itr = g->begin();
    agi::GraphVertex* vert;
    while ((vert = g->iterate(itr)));
    
    agi::destroyGraph(g);
  }
  else {
    agi::Ngraph* g = buildHyperGraph();

    agi::EdgeIterator* itr = g->begin(0);
    agi::GraphEdge* edge;
    while ((edge = g->iterate(itr)));
    g->destroy(itr);  

    agi::destroyGraph(g);

  }
  PCU_Switch_Comm(MPI_COMM_WORLD);
  PCU_Barrier();
  if (!old_self) 
    printf("All Tests Passed\n");

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}
