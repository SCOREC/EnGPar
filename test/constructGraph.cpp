#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include <set>
#include "buildGraphs.h"

void testGraph();
void testHyperGraph();
void testGraphParts();
void testHyperGraphParts();
void testRequirements();
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  testGraph();

  PCU_Barrier();
  
  testHyperGraph();

  PCU_Barrier();

  testGraphParts();

  PCU_Barrier();

  testHyperGraphParts();

  PCU_Barrier();

  testRequirements();

  PCU_Barrier();

  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n");
  
  EnGPar_Finalize();
  MPI_Finalize();
}

//Tests a ring of vertices
//  where each part gets 4 continuous vertices
void testGraph() {
  if (!PCU_Comm_Self())
    printf("Testing Regular Graph\n");
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
  int total =0;
  while ((gv = graph->iterate(g_itr))) {
    assert(graph->localID(gv)>=graph->numLocalVtxs());
    printf("%d %ld %ld %d\n",PCU_Comm_Self(),graph->localID(gv),graph->globalID(gv),graph->owner(gv));
    assert(graph->owner(gv)!=PCU_Comm_Self());
    total++;
  }
  assert(total==graph->numGhostVtxs());
  
  agi::destroyGraph(graph);

}

void testHyperGraph() {
  if (!PCU_Comm_Self())
    printf("Testing HyperGraph\n");
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
    agi::GraphEdge* e;
    while ((e = graph->iterate(eitr))) {
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
}


//Tests a ring of vertices
//  where each part gets 4 continuous vertices
void testGraphParts() {
  if (!PCU_Comm_Self())
    printf("Testing Regular Graph Parts\n");
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
}


void testHyperGraphParts() {
  if (!PCU_Comm_Self())
    printf("Testing HyperGraph Parts\n");
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
    agi::GraphEdge* e;
    while ((e = graph->iterate(eitr))) {
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
}


void testRequirements() {
  if (!PCU_Comm_Self())
    printf("Testing Requirements Graph\n");  
  agi::Ngraph* g = buildRequirementsGraph();
  agi::checkValidity(g);
  agi::destroyGraph(g);
}
