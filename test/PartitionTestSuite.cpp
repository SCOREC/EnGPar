#include <ngraph.h>
#include <engpar_support.h>
#include <iostream>
#include <PCU.h>
#include <vector>
#include "buildGraphs.h"
#include "gatherGraphs.h"
#include "TestingSuite.h"
#include <engpar.h>

int testVtxBalancer(agi::Ngraph*);
int testBalancer(agi::Ngraph*);
int testMultipleBalances(agi::Ngraph*);

int testWeightBalancer_1();
int testWeightBalancer_4();
int testWeightBalancer_100();

int testGlobalSplit(agi::Ngraph*);
int testLocalSplit(agi::Ngraph*);
int testSplitAndBalance(agi::Ngraph*);

void switchToOriginals(int smallSize, bool& isOriginal, MPI_Comm& newComm);

#define SPLIT_TEST 4

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if (argc==1) {
    if (!PCU_Comm_Self())
      EnGPar_Warning_Message("Usage: %s <operation> [trial number]\n"
                             "    operation 0 = vtxBalance\n"
                             "    operation 1 = balance\n"
                             "    operation 2 = balanceMultiple\n"
                             "    operation 3 = balanceWeights\n"
                             "    operation %d = ParMETIS global split\n"
                             "    operation %d = ParMETIS local split\n"
                             "    operation %d = ParMETIS split and balance\n"
                             ,argv[0],SPLIT_TEST,SPLIT_TEST+1,SPLIT_TEST+2);

    EnGPar_Finalize();
    MPI_Finalize();
    return 1;
  }
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
  if (operation==3) {
    suite.addFineTest("Weight Balancer with 1 Component",testWeightBalancer_1);
    suite.addFineTest("Weight Balancer with 4 Components",testWeightBalancer_4);
    suite.addFineTest("Weight Balancer with 100 Components",testWeightBalancer_100);
  }

  bool isOriginal = true;
  if (operation>=4) {
    MPI_Comm newComm;
    switchToOriginals(2, isOriginal,newComm);
    //Switch the internal communicator (this changes PCU so use PCU_Comm_... with caution)
    EnGPar_Switch_Comm(newComm);
  }

  //Gather the graphs for the general tests
  if (isOriginal) {
    gatherEBINGraphs(suite);
    gatherBGDGraphs(suite);
  }
  EnGPar_Switch_Comm(MPI_COMM_WORLD);
  suite.fillEmptyTestGraphs();
  
  //Gather general tests that run on the graphs collected prior to this
  if (operation==0)
    suite.addGeneralTest("Vertex Balancer",testVtxBalancer);
  else if (operation==1)
    suite.addGeneralTest("General Balancer",testBalancer);
  else if (operation==2)
    suite.addGeneralTest("Balance Twice",testMultipleBalances);
  else if (operation == SPLIT_TEST)
    suite.addGeneralTest("Global ParMETIS split",testGlobalSplit);
  else if (operation == SPLIT_TEST+1)
    suite.addGeneralTest("Local ParMETIS split",testLocalSplit);
  else if (operation == SPLIT_TEST+2)
    suite.addGeneralTest("Global Split and MC Balance",testSplitAndBalance);
  
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
  engpar::balanceWeights(input,-1);
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
  engpar::balanceWeights(input,-1);
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
  engpar::balanceWeights(input,-1);
  //TODO: add balancing checks
  agi::checkValidity(g);
  agi::destroyGraph(g);

  return 0;
}

int testVtxBalancer(agi::Ngraph* g) {

  //Create the balancer
  engpar::balanceVertices(g, 1.1, .1, -1);

  agi::checkValidity(g);
  
  if (engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,-1)) >= 1.11)
    return 1;
  return 0;

}
int testBalancer(agi::Ngraph* g) {
  double step_factor = 0.1;
  engpar::DiffusiveInput* input = engpar::createDiffusiveInput(g,step_factor);
  for (agi::etype t = 0; t < g->numEdgeTypes(); t++)
    input->addPriority(t,1.1);
  input->addPriority(-1,1.1);

  input->maxIterationsPerType=50;
  input->maxIterations=75;

  //Create the balancer
  engpar::balance(input,-1);

  //Ensure the graph is still valid
  agi::checkValidity(g);

  //Ensure the graph was balanced to the target tolerance
  for (agi::etype t = 0; t < g->numEdgeTypes(); t++)
    if (engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,t)) >= 1.11)
      return t+2;
  if (engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,-1)) >= 1.11)
    return 1;
  return 0;
}
int testMultipleBalances(agi::Ngraph* g) {

  double step_factor = 0.1;
  engpar::Input* input = engpar::createDiffusiveInput(g,step_factor);
  input->addPriority(-1,1.1);

  //Create the balancer
  engpar::balance(input,-1);

  //Ensure the graph is still valid
  agi::checkValidity(g);

  if (engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,-1)) >= 1.11)
    return 1;
  
  //Alter weights so odd processes are heavier than even processes
  agi::wgt_t* weights = new agi::wgt_t[g->numLocalVtxs()];
  agi::wgt_t w = PCU_Comm_Self()%2+1;

  for (agi::lid_t i =0;i<g->numLocalVtxs();i++) {
    weights[i] = w;
  }
  g->setWeights(weights);
  delete [] weights;

  g->resetOwnership();
  
  input = engpar::createDiffusiveInput(g,step_factor);
  input->addPriority(-1,1.1);

  //Create the balancer
  engpar::balance(input,-1);

  //Ensure the graph is still valid
  agi::checkValidity(g);

  if (engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,-1)) >= 1.11)
    return 2;
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

int testGlobalSplit(agi::Ngraph* g) {
  //Application code:
  bool isOriginal = g->numLocalVtxs()>0;
  MPI_Comm newComm;
  int split_factor = 2;
  switchToOriginals(split_factor, isOriginal, newComm);

  //Switch the internal communicator (this changes PCU so use PCU_Comm_... with caution)
  EnGPar_Switch_Comm(newComm);

  //Create the input
  double tolerance = 1.05;
  agi::etype t = 0;
  engpar::Input* input = engpar::createGlobalSplitInput(g,newComm,MPI_COMM_WORLD, isOriginal,
                                                        tolerance,t);

  engpar::split(input,engpar::GLOBAL_PARMETIS);

  if (PCU_Get_Comm() != MPI_COMM_WORLD)
    return 2;

  if (g->numGlobalVtxs() / PCU_Comm_Peers() >= 20 && engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,-1)) >= 1.11)
    return 1;
  
  //Application continues:
  MPI_Comm_free(&newComm);
  return 0;
}

int testLocalSplit(agi::Ngraph* g) {
  
  //Application code:
  int split_factor = 2;
  bool isOriginal = g->numLocalVtxs()>0;
  MPI_Comm newComm;
  switchToOriginals(split_factor, isOriginal, newComm);

  //Switch the internal communicator (this changes PCU so use PCU_Comm_... with caution)
  EnGPar_Switch_Comm(newComm);

  //Create the input

  double tolerance = 1.05;
  agi::etype t = 0;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  agi::part_t* others = new agi::part_t[split_factor];
  for (int i=0;i<split_factor;i++) {
    others[i] = my_rank+i;
  }

  engpar::Input* input = engpar::createLocalSplitInput(g,newComm,MPI_COMM_WORLD, isOriginal,
                                                       2,tolerance,others,t);

  engpar::split(input,engpar::LOCAL_PARMETIS);

  if (PCU_Get_Comm() != MPI_COMM_WORLD)
    return 2;

  if (g->numGlobalVtxs() / PCU_Comm_Peers() >= 20 && engpar::EnGPar_Get_Imbalance(engpar::getWeight(g,-1)) >= 1.11)
    return 1;
  
  //Application continues:
  MPI_Comm_free(&newComm);
  
  return 0;
}
int testSplitAndBalance(agi::Ngraph* g) {
  int ierr = testGlobalSplit(g);
  if (ierr>0)
    return ierr;
  ierr = testBalancer(g);
  if (ierr>0)
    return ierr+10;
  
  return 0;
}
