#include <ngraph.h>
#include <engpar_support.h>
#include <iostream>
#include <PCU.h>
#include <vector>
#include "buildGraphs.h"
#include "gatherGraphs.h"
#include "TestingSuite.h"
#include <engpar.h>

int testWeightBalancer_1();
int testWeightBalancer_4();
int testWeightBalancer_100();


int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  int operation = -1;
  if (argc>1)
    operation = atoi(argv[1]);
  int trial = -1;
  if (argc>2)
    trial = atoi(argv[2]);

  
  EnGPar_Set_Verbosity(-1);
  //Create the testing suite
  TestingSuite suite(argv[0]);
    
  //Gather specific tests that have more fine grain checks
  //suite.addFineTest("Switch Comm", switchComm);
  if (operation==6) {
    suite.addFineTest("WeightBalancer with 1 Component",testWeightBalancer_1);
    suite.addFineTest("WeightBalancer with 4 Components",testWeightBalancer_4);
    suite.addFineTest("WeightBalancer with 100 Components",testWeightBalancer_100);
  }
  
  //Gather the graphs for the general tests
  gatherBuildGraphs(suite);
  gatherEBINGraphs(suite);
  gatherBGDGraphs(suite);
  
  //Gather general tests that run on the graphs collected prior to this
  //suite.addGeneralTest("Traverse Edges",traverseEdges);
  
  //Run the tests and get the number of failures
  int ierr = suite.runTests(trial);

  suite.deleteTestGraphs();
  EnGPar_Finalize();
  MPI_Finalize();
  return ierr;
}

agi::Ngraph* makeArtificialXGCMGraph(int numComps) {
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
  delete [] verts;
  delete [] weights;
  
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
  delete [] edge_ids;
  delete [] degs;
  delete [] pins_to_verts;
  
  std::unordered_map<agi::gid_t,agi::part_t> owns;
  for (int i=0;i<num_edges;i++) {
    for (int j=0;j<PCU_Comm_Peers();j++) {
      if (j!=PCU_Comm_Self())
        owns[i+j*num_edges] = j;
    }
  }
  g->constructGhosts(owns);
  return g;
}

int testWeightBalancer_1() {
  double tol = 1.025;
  agi::Ngraph* g = makeArtificialXGCMGraph(1);
  double step_factor = .25;
  engpar::WeightInput* input = engpar::createWeightInput(g,tol,step_factor,0);
  engpar::balanceWeights(input,1);
  //TODO: add balancing checks
  agi::checkValidity(g);
  agi::destroyGraph(g);
  return 0;
}
int testWeightBalancer_4() {
  double tol = 1.025;
  agi::Ngraph* g = makeArtificialXGCMGraph(4);
  double step_factor = .25;
  engpar::WeightInput* input = engpar::createWeightInput(g,tol,step_factor,0);
  engpar::balanceWeights(input,1);
  //TODO: add balancing checks
  agi::checkValidity(g);
  agi::destroyGraph(g);
  return 0;
}
int testWeightBalancer_100() {
  double tol = 1.025;
  agi::Ngraph* g = makeArtificialXGCMGraph(100);
  double step_factor = .25;
  engpar::WeightInput* input = engpar::createWeightInput(g,tol,step_factor,0);
  engpar::balanceWeights(input,1);
  //TODO: add balancing checks
  agi::checkValidity(g);
  agi::destroyGraph(g);

  return 0;
}
