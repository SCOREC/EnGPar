#include <ngraph.h>
#include <engpar_support.h>
#include <binGraph.h>
#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include <PCU.h>
#include "buildGraphs.h"
#include "gatherGraphs.h"
#include "TestingSuite.h"

int testGraph();
int testHyperGraph();
int testGraphParts();
int testHyperGraphParts();
int testBuildGraphs();
int testEBINGraphs();
int testBGDGraphs();

int isValid(agi::Ngraph*);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Set_Verbosity(-1);

  double mem = EnGPar_Peak_Memory();
  if (!PCU_Comm_Self())
    printf("Peak memory before tests is %.4f\n", mem);
  //Create the testing suite
  TestingSuite suite(argv[0]);

  mem = EnGPar_Peak_Memory();
  if (!PCU_Comm_Self())
    printf("Peak memory after creating suite is %.4f\n", mem);

  //Gather specific tests that have more fine grain checks
  suite.addFineTest("Construct Graph",testGraph);
  suite.addFineTest("Construct HyperGraph",testHyperGraph);
  suite.addFineTest("Construct Graph in Parts",testGraphParts);
  suite.addFineTest("Construct HyperGraph in Parts",testHyperGraphParts);
  suite.addFineTest("Construct Build Graphs",testBuildGraphs);
  #ifndef ENGPAR_BIG_ENDIAN
  suite.addFineTest("Construct .ebin Graphs",testEBINGraphs);
  #endif
  suite.addFineTest("Construct .bgd Graphs",testBGDGraphs);


  mem = EnGPar_Peak_Memory();
  if (!PCU_Comm_Self())
    printf("Peak memory after creating tests is %.4f\n", mem);

  //Gather the graphs for the general tests
  gatherBuildGraphs(suite);
  gatherEBINGraphs(suite);
  gatherBGDGraphs(suite);

  mem = EnGPar_Peak_Memory();
  if (!PCU_Comm_Self())
    printf("Peak memory after creating graphs is %.4f\n", mem);

  //Gather general tests that run on the graphs collected prior to this
  suite.addGeneralTest("Check Validity",isValid);
  
  //Run the tests and get the number of failures
  int ierr = suite.runTests();

  mem = EnGPar_Peak_Memory();
  if (!PCU_Comm_Self())
    printf("Peak memory after running tests is %.4f\n", mem);
  
  suite.deleteTestGraphs();
  EnGPar_Finalize();
  MPI_Finalize();
  return ierr;
}


//Tests a ring of vertices
//  where each part gets 4 continuous vertices
int testGraph() {
  agi::Ngraph* graph = buildGraph();
  agi::lid_t local_verts = 4;
  agi::gid_t global_verts = 4*PCU_Comm_Peers();

  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==2);
  }
  else
    assert(graph->numGhostVtxs()==0);
  assert(graph->numLocalEdges()==local_verts*2);
  assert(graph->numEdgeTypes()==1);
  assert(!graph->isHyper());

  std::set<agi::gid_t> vs;
  for (int i = -1 ; i<(int)local_verts+1;i++)  {
    agi::gid_t vert = (PCU_Comm_Self()*local_verts+global_verts+i)%global_verts;
    vs.insert(vert);
  }
  agi::GraphVertex* v;
  agi::VertexIterator* vitr = graph->begin();
  while ((v = graph->iterate(vitr))) {
    assert(vs.find(graph->globalID(v))!=vs.end());
    assert(graph->localID(v)<graph->numTotalVtxs());
    
    if (graph->localID(v)>graph->numLocalVtxs())
      assert(graph->owner(v)!=PCU_Comm_Self());
    agi::GraphVertex* other;
    agi::GraphIterator* gitr = graph->adjacent(v);
    while ((other = graph->iterate(gitr))) {
      assert(vs.find(graph->globalID(other))!=vs.end());
      agi::GraphEdge* edge = graph->edge(gitr);
      assert(graph->degree(edge)==2);
      agi::PinIterator* pitr = graph->pins(edge);
      agi::GraphVertex* v1 = graph->u(edge);
      agi::GraphVertex* v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(v)==graph->localID(v2));
      v1 = graph->v(edge);
      v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(other)==graph->localID(v2));

      assert(!graph->iterate(pitr));
      graph->destroy(pitr);
    }
    graph->destroy(gitr);
  }

  agi::GhostIterator* g_itr = graph->beginGhosts();
  agi::GraphVertex* gv;
  int total = 0;
  while ((gv = graph->iterate(g_itr))) {
    assert(graph->localID(gv)>=graph->numLocalVtxs());
    assert(graph->owner(gv)!=PCU_Comm_Self());
    total++;
  }
  assert(total==graph->numGhostVtxs());
  
  agi::destroyGraph(graph);
  return 0;
}

int testHyperGraph() {
  agi::lid_t local_verts = 4;
  agi::lid_t local_edges = 3;

  agi::Ngraph* graph = buildHyperGraph();
  
  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==4);
  }
  assert(graph->numLocalEdges()==local_edges+(PCU_Comm_Peers()>1));
  assert(graph->numEdgeTypes()==1);
  assert(graph->isHyper());
  //assert(graph->numLocalPins()==pins.size());

  agi::lid_t vert_degs[4] = {2,3,2,3};
  agi::lid_t edge_degs[4] = {4,2,4,4};
  agi::lid_t ghost_degs[4] = {0,0,2,2};
  agi::gid_t pin[14] = {0,1,2,3,1,3,2,3,4,5,6,7,0,1};
  agi::VertexIterator* itr = graph->begin();
  agi::GraphVertex* vtx;
  agi::lid_t i=0;
  while ((vtx = graph->iterate(itr))) {
    assert(graph->localID(vtx)<graph->numLocalVtxs());
    assert(graph->globalID(vtx)<graph->numGlobalVtxs());
    assert(graph->degree(vtx)==vert_degs[i]);
    agi::lid_t count = 0;
    agi::EdgeIterator* eitr = graph->edges(vtx);
    while ((graph->iterate(eitr))) {
      count++;
    }
    assert(count==vert_degs[i]);
    graph->destroy(eitr);
    i++;
  }
  assert(i==graph->numLocalVtxs());
  i=0;
  
  agi::GraphEdge* e;
  agi::EdgeIterator* eitr = graph->begin(0);
  int k=0;
  while ((e = graph->iterate(eitr))) {
    agi::lid_t deg = graph->degree(e);
    assert(deg==edge_degs[i]);
    agi::lid_t count =0;
    agi::lid_t ghost =0;
    agi::PinIterator* pitr = graph->pins(e);
    for (agi::lid_t j=0;j<deg;j++) {
      vtx = graph->iterate(pitr);
      if (PCU_Comm_Peers()>1||pin[k]<graph->numLocalVtxs())
        assert(graph->localID(vtx)==pin[k++]);
      assert(graph->localID(vtx)<graph->numTotalVtxs());
      assert(graph->globalID(vtx)<graph->numGlobalVtxs());
      count++;
      if (graph->localID(vtx)>=graph->numLocalVtxs()) {
        assert(PCU_Comm_Peers()>1);
        assert(graph->owner(vtx)!=PCU_Comm_Self());
        ghost++;
      }
    }
    graph->destroy(pitr);
    assert(ghost==ghost_degs[i]*(PCU_Comm_Peers()>1));
    assert(count==edge_degs[i]);
    i++;
  }
  graph->destroy(eitr);

  agi::GhostIterator* g_itr = graph->beginGhosts();
  agi::GraphVertex* gv;
  int total =0;
  while ((gv = graph->iterate(g_itr))) {
    assert(graph->localID(gv)>=graph->numLocalVtxs());
    assert(graph->owner(gv)!=PCU_Comm_Self());
    total++;
  }
  assert(total==graph->numGhostVtxs());

  agi::destroyGraph(graph);
  return 0;
}


//Tests a ring of vertices
//  where each part gets 4 continuous vertices
int testGraphParts() {
  agi::Ngraph* graph  = buildGraphParts();
  agi::lid_t local_verts = 4;
  agi::gid_t global_verts = 4*PCU_Comm_Peers();
  agi::etype t = 0;
  agi::etype t2 = 1;
  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==2+(PCU_Comm_Self()!=0));
  }
  else
    assert(graph->numGhostVtxs()==0);
  assert(graph->numLocalEdges(t)==local_verts*2);
  assert(graph->numLocalEdges(t2) == local_verts+(PCU_Comm_Self()!=0));
  assert(graph->numEdgeTypes()==2);
  assert(!graph->isHyper());

  std::set<agi::gid_t> vs;
  for (int i = -1 ; i<(int)local_verts+1;i++)  {
    agi::gid_t vert = (PCU_Comm_Self()*local_verts+global_verts+i)%global_verts;
    vs.insert(vert);
  }
  if (PCU_Comm_Self())
    vs.insert(1);
  agi::GraphVertex* v;
  agi::VertexIterator* vitr = graph->begin();
  while ((v = graph->iterate(vitr))) {
    assert(vs.find(graph->globalID(v))!=vs.end());
    assert(graph->localID(v)<graph->numTotalVtxs());
    
    if (graph->localID(v)>graph->numLocalVtxs())
      assert(graph->owner(v)!=PCU_Comm_Self());
    agi::GraphVertex* other;
    agi::GraphIterator* gitr = graph->adjacent(v,t);
    while ((other = graph->iterate(gitr))) {
      assert(vs.find(graph->globalID(other))!=vs.end());
      agi::GraphEdge* edge = graph->edge(gitr);
      assert(graph->degree(edge)==2);
      agi::PinIterator* pitr = graph->pins(edge);
      agi::GraphVertex* v1 = graph->u(edge);
      agi::GraphVertex* v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(v)==graph->localID(v2));
      v1 = graph->v(edge);
      v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(other)==graph->localID(v2));

      assert(!graph->iterate(pitr));
      graph->destroy(pitr);
    }
    graph->destroy(gitr);
    gitr = graph->adjacent(v,t2);
    while ((other = graph->iterate(gitr))) {
      assert(vs.find(graph->globalID(other))!=vs.end());
      agi::GraphEdge* edge = graph->edge(gitr);
      assert(graph->degree(edge)==2);
      agi::PinIterator* pitr = graph->pins(edge);
      agi::GraphVertex* v1 = graph->u(edge);
      agi::GraphVertex* v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(v)==graph->localID(v2));
      v1 = graph->v(edge);
      v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(other)==graph->localID(v2));
      
      assert(!graph->iterate(pitr));
      graph->destroy(pitr);
    }
    graph->destroy(gitr);
  }

  //Remove the edges of type t and then check if t2 still works
  graph->removeEdges(t);
  assert(graph->numLocalEdges(t)==0);
  assert(graph->numGlobalEdges(t)==0);
  vitr = graph->begin();
  while ((v = graph->iterate(vitr))) {
    agi::GraphVertex* other;
    agi::GraphIterator* gitr = graph->adjacent(v,t2);
    while ((other = graph->iterate(gitr))) {
      assert(vs.find(graph->globalID(other))!=vs.end());
      agi::GraphEdge* edge = graph->edge(gitr);
      assert(graph->degree(edge)==2);
      agi::PinIterator* pitr = graph->pins(edge);
      agi::GraphVertex* v1 = graph->u(edge);
      agi::GraphVertex* v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(v)==graph->localID(v2));
      v1 = graph->v(edge);
      v2 = graph->iterate(pitr);
      assert(graph->localID(v1)==graph->localID(v2));
      assert(graph->localID(other)==graph->localID(v2));
      
      assert(!graph->iterate(pitr));
      graph->destroy(pitr);
    }
    graph->destroy(gitr);
  }
  
  agi::destroyGraph(graph);
  return 0;
}


int testHyperGraphParts() {
  agi::Ngraph* graph = buildHyperGraphParts();
  agi::lid_t local_verts = 4;
  agi::lid_t local_edges = 3;

  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==4);
  }
  assert(graph->numLocalEdges()==local_edges+(PCU_Comm_Peers()>1));
  assert(graph->numEdgeTypes()==1);
  assert(graph->isHyper());
  //assert(graph->numLocalPins()==pins.size());

  agi::lid_t vert_degs[4] = {2,3,2,3};
  agi::lid_t edge_degs[4] = {4,2,4,4};
  agi::lid_t ghost_degs[4] = {0,0,2,2};
  agi::gid_t pin[14] = {0,1,2,3,1,3,2,3,4,5,6,7,0,1};
  agi::VertexIterator* itr = graph->begin();
  agi::GraphVertex* vtx;
  agi::lid_t i=0;
  while ((vtx = graph->iterate(itr))) {
    assert(graph->localID(vtx)<graph->numLocalVtxs());
    assert(graph->globalID(vtx)<graph->numGlobalVtxs());
    assert(graph->degree(vtx)==vert_degs[i]);
    agi::lid_t count = 0;
    agi::EdgeIterator* eitr = graph->edges(vtx);
    while ((graph->iterate(eitr))) {
      count++;
    }
    assert(count==vert_degs[i]);
    graph->destroy(eitr);
    i++;
  }
  assert(i==graph->numLocalVtxs());
  i=0;
  
  agi::GraphEdge* e;
  agi::EdgeIterator* eitr = graph->begin(0);
  int k=0;
  while ((e = graph->iterate(eitr))) {
    agi::lid_t deg = graph->degree(e);
    assert(deg==edge_degs[i]);
    agi::lid_t count =0;
    agi::lid_t ghost =0;
    agi::PinIterator* pitr = graph->pins(e);
    for (agi::lid_t j=0;j<deg;j++) {
      vtx = graph->iterate(pitr);
      if (PCU_Comm_Peers()>1||pin[k]<graph->numLocalVtxs())
        assert(graph->localID(vtx)==pin[k++]);
      assert(graph->localID(vtx)<graph->numTotalVtxs());
      assert(graph->globalID(vtx)<graph->numGlobalVtxs());
      count++;
      if (graph->localID(vtx)>=graph->numLocalVtxs()) {
        assert(PCU_Comm_Peers()>1);
        assert(graph->owner(vtx)!=PCU_Comm_Self());
        ghost++;
      }
    }
    graph->destroy(pitr);
    assert(ghost==ghost_degs[i]*(PCU_Comm_Peers()>1));
    assert(count==edge_degs[i]);
    i++;
  }
  graph->destroy(eitr);
  agi::destroyGraph(graph);
  return 0;
}

int testBuildGraphs() {
  agi::Ngraph* g = buildGraph();
  agi::destroyGraph(g);
  g = buildHyperGraph();
  agi::destroyGraph(g);
  g = buildGraphParts();
  agi::destroyGraph(g);
  g = buildHyperGraphParts();
  agi::destroyGraph(g);
  g = buildRequirementsGraph();
  agi::destroyGraph(g);
  g = buildEmptyGraph();
  agi::destroyGraph(g);
  if (PCU_Comm_Peers()==2) {
    g = buildDisconnected2Graph();
    agi::destroyGraph(g);
    g = buildUnbalancedLineHG();
    agi::destroyGraph(g);
  }
  return 0;

}
int testEBINGraphs() {
  char file[1024];
  //Ring graph
  sprintf(file,"%s/ring.ebin",ENGPAR_GRAPHS);
  agi::Ngraph* g = agi::createBinGraph(file);
  agi::destroyGraph(g);

  //Tree graph
  sprintf(file,"%s/tree.ebin",ENGPAR_GRAPHS);
  g = agi::createBinGraph(file);
  agi::destroyGraph(g);

  //Gnutella graph
  sprintf(file,"%s/gnutella.ebin",ENGPAR_GRAPHS);
  g = agi::createBinGraph(file);
  agi::destroyGraph(g);
  return 0;
}
int testBGDGraphs() {
  //Cube tests (comm size = 1 or 4)
  char file[1024];
  file[0] = '\0';
  if (PCU_Comm_Peers() == 1) {
    sprintf(file,"%s/cube/cube",ENGPAR_GRAPHS);
  }
  else if (PCU_Comm_Peers()==4) {
    sprintf(file,"%s/cube/4/",ENGPAR_GRAPHS);
  }
  if (file[0]!='\0') {
    agi::Ngraph* g = agi::createEmptyGraph();
    g->loadFromFile(file);
    agi::destroyGraph(g);
  }

  if (PCU_Comm_Peers() == 4) {
    sprintf(file,"%s/torus/4/",ENGPAR_GRAPHS);
    agi::Ngraph* g = agi::createEmptyGraph();
    g->loadFromFile(file);
    agi::destroyGraph(g);
    sprintf(file,"%s/torus/4_01/",ENGPAR_GRAPHS);
    g = agi::createEmptyGraph();
    g->loadFromFile(file);
    agi::destroyGraph(g);
  }
  return 0;
}

int isValid(agi::Ngraph* g) {
  return (!agi::checkValidity(g));
}
