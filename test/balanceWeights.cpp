#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <cstring>
#include "buildGraphs.h"


int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  
  if (argc!= 3) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <numberOfComponents> <tolerance>\n", argv[0]);
    }
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }
  int numComps = atoi(argv[1]);
  double tol = atof(argv[2]);
  
  //Create the graph where each process has 1 vertex per component
  // each vertex has weight = PCU_Comm_Self()*100
  /* The Ngraph for parts=4, numComps = 3 would be:
       numbers=verts, letters = hyperedges
       
       1   4   2   5   3   6
        \ /     \ /     \ / 
         A       B       C
        / \     / \     / \
       7   10  8   11  9   12
  */
  agi::Ngraph* g= buildEmptyGraph();
  agi::lid_t num_verts = numComps;
  agi::gid_t* verts = new agi::gid_t[num_verts];
  agi::wgt_t* weights = new agi::wgt_t[num_verts];
  for (int i=0;i<num_verts;i++) {
    verts[i] = PCU_Comm_Self()*numComps+i;
    weights[i] = PCU_Comm_Self()*100;
  }
  g->constructVerts(true,num_verts,verts,weights);

  agi::gid_t num_edges = numComps;
  agi::gid_t* edge_ids = new agi::gid_t[num_edges];
  agi::lid_t* degs = new agi::lid_t[num_edges];
  agi::gid_t* pins_to_verts = new agi::gid_t[num_edges*PCU_Comm_Peers()];
  for (int i=0;i<num_edges;i++) {
    edge_ids[i] = i;
    degs[i] = PCU_Comm_Peers();
    for (int j=0;j<PCU_Comm_Peers();j++) {
      pins_to_verts[i*PCU_Comm_Peers()+j] = i+j*num_edges;
    }
  }
  g->constructEdges(num_edges,edge_ids,degs,pins_to_verts);

  std::unordered_map<agi::gid_t,agi::part_t> owns;
  for (int i=0;i<num_edges;i++) {
    for (int j=0;j<PCU_Comm_Peers();j++) {
      if (j!=PCU_Comm_Self())
        owns[i+j*num_edges] = j;
    }
  }
  g->constructGhosts(owns);
  
  engpar::evaluatePartition(g);

  double step_factor = 0.1;
  engpar::WeightInput* input = engpar::createWeightInput(g,1.05,step_factor,0);

  //Create the balancer
  engpar::balanceWeights(input,0);

  //Evaluate the new partition
  engpar::evaluatePartition(g);

  //Ensure the graph is still valid
  agi::checkValidity(g);
  
  //Destroy graph
  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}
