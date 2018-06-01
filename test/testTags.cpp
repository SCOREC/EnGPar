#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include <set>
#include "buildGraphs.h"
void tagGraph();

void tagHyperGraph();

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  tagGraph();
  
  PCU_Barrier();

  tagHyperGraph();
  
  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n");

  EnGPar_Finalize();
  MPI_Finalize();

  
}

void tagGraph() {
  if (!PCU_Comm_Self())
    printf("Tagging regular graph\n");
  agi::Ngraph* g = buildGraphParts();

  agi::GraphTag* tag = g->createIntTag();
  agi::VertexIterator* itr= g->begin();
  agi::GraphVertex* vtx;

  int i=0;
  while ((vtx = g->iterate(itr))) {
    g->setIntTag(tag,vtx,i);
    i++;
  }

  itr = g->begin();
  i=0;
  while ((vtx = g->iterate(itr))) {
    assert(g->getIntTag(tag,vtx)==i);
    i++;
  }

  g->destroyTag(tag);

  tag = g->createLongTag(0);
  agi::GraphEdge* edge;
  agi::EdgeIterator* eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    g->setDoubleTag(tag,edge,g->degree(edge));
  }
  g->destroy(eitr);
  eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    assert(g->getDoubleTag(tag,edge)==2);
  }
  g->destroy(eitr);
  g->destroyTag(tag);

  tag = g->createLongTag(1);

  eitr = g->begin(1);
  while ((edge=g->iterate(eitr))) {
    g->setDoubleTag(tag,edge,g->degree(edge));
  }
  g->destroy(eitr);
  eitr = g->begin(1);
  while ((edge=g->iterate(eitr))) {
    assert(g->getDoubleTag(tag,edge)==2);
  }
  g->destroy(eitr);

  agi::destroyGraph(g);
  
}

void tagHyperGraph() {
  if (!PCU_Comm_Self())
    printf("Tagging regular graph\n");

  agi::Ngraph* g = buildHyperGraphParts();

  agi::GraphTag* tag = g->createIntTag();
  agi::VertexIterator* itr= g->begin();
  agi::GraphVertex* vtx;

  int i=0;
  while ((vtx = g->iterate(itr))) {
    g->setIntTag(tag,vtx,i);
    i++;
  }

  itr = g->begin();
  i=0;
  while ((vtx = g->iterate(itr))) {
    assert(g->getIntTag(tag,vtx)==i);
    i++;
  }

  g->destroyTag(tag);

  tag = g->createLongTag(0);
  agi::GraphEdge* edge;
  agi::EdgeIterator* eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    g->setDoubleTag(tag,edge,g->degree(edge));
  }
  g->destroy(eitr);
  eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    assert(g->getDoubleTag(tag,edge)==g->degree(edge));
  }
  g->destroy(eitr);
  g->destroyTag(tag);

  agi::destroyGraph(g);
}
