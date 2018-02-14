#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include <set>
#include "buildGraphs.h"

void testGraphParts();
void testHyperGraphParts();
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  assert(PCU_Comm_Peers()>1);
  
  testGraphParts();

  PCU_Barrier();

  testHyperGraphParts();

  PCU_Barrier();


  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n");
  
  EnGPar_Finalize();
  MPI_Finalize();
}


//Tests a ring of vertices
//  where each part gets 4 continuous vertices
void testGraphParts() {
  if (!PCU_Comm_Self())
    printf("Testing Regular Graph Parts\n");
  agi::Ngraph* graph  = buildGraphParts();

  agi::part_t partition[4];

  for (int i=0;i<2;i++)
    partition[i] = PCU_Comm_Self();
  
  for (int i=2;i<4;i++)
    partition[i] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();

  graph->repartition(partition);

  checkValidity(graph);

  agi::destroyGraph(graph);
}


void testHyperGraphParts() {
  if (!PCU_Comm_Self())
    printf("Testing HyperGraph Parts\n");
  agi::Ngraph* graph = buildHyperGraphParts();

  agi::part_t partition[4];

  for (int i=0;i<2;i++)
    partition[i] = PCU_Comm_Self();
  
  for (int i=2;i<4;i++)
    partition[i] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();

  graph->repartition(partition);

  checkValidity(graph);


  agi::destroyGraph(graph);
}
