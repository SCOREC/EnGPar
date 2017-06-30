#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include <set>
#include "buildGraphs.h"

void migrateGraphParts();

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  agi::Ngraph* g = buildGraph();
  PCU_Barrier();
  agi::Migration* plan = new agi::Migration;
  agi::GraphVertex* v;
  agi::VertexIterator* itr=g->begin();
  while ((v = g->iterate(itr))) {
    agi::GraphIterator* gitr = g->adjacent(v);
    agi::GraphVertex* other;
    while ((other=g->iterate(gitr))) {
      int owner = g->owner(other);
      if (owner>PCU_Comm_Self()) {
        plan->insert(std::make_pair(v,owner));
        break;
      }
    }
  }
  g->setOriginalOwners();
  if (!PCU_Comm_Self())
    printf("Migrating traditional graph\n");
  g->migrate(plan);
  
  agi::checkValidity(g);
  agi::destroyGraph(g);
  
  PCU_Barrier();

  g = buildHyperGraph();
  PCU_Barrier();
  plan = new agi::Migration;

  itr=g->begin();
  v = g->iterate(itr);
  (*plan)[v] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();
  g->setOriginalOwners();
  if (!PCU_Comm_Self())
    printf("Migrating hypergraph\n");
  g->migrate(plan);
  agi::checkValidity(g);
  agi::destroyGraph(g);
  PCU_Barrier();

  migrateGraphParts();

  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n");
  

  EnGPar_Finalize();
  MPI_Finalize();

}
void migrateGraphParts() {
  if (!PCU_Comm_Self())
    printf("Testing Migrating Regular Graph Parts\n");
  agi::Ngraph* graph = buildGraphParts();
  agi::etype t = 0;
  agi::etype t2 = 1;
  agi::lid_t local_verts = 4;
  agi::gid_t global_verts = 4*PCU_Comm_Peers();

  //Create a migration and send nothing to test if the reconstructed graph
  // is the same even with 2 edge types
  graph->setOriginalOwners();
  agi::Migration* plan = new agi::Migration;
  graph->migrate(plan);  
  
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
  agi::destroyGraph(graph);
}
