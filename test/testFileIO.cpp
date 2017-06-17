#include <apfGraph.h>
#include <apfMesh2.h>
#include <cassert>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <cstdlib>
#include <stdint.h>
#include <apfNumbering.h>
#include <cstring>
#include <engpar_support.h>

void testSizes(agi::Ngraph*,agi::Ngraph*);
void testVertices(agi::Ngraph*,agi::Ngraph*);
void testEdges(agi::Ngraph*,agi::Ngraph*);

void testGraphs(agi::Ngraph* g1,agi::Ngraph* g2) {
  
  testSizes(g1,g2);
  testVertices(g1,g2);  
  testEdges(g1,g2);  
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  EnGPar_Initialize();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <save prefix>\n", argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  agi::Ngraph* g = agi::createAPFGraph(m,3,2);
  std::vector<agi::wgt_t> wgts;
  g->setEdgeWeights(wgts,0);
  PCU_Barrier();

  g->saveToFile(argv[3]);
  agi::Ngraph* gLoad = new agi::Ngraph;
  gLoad->loadFromFile(argv[3]);
  
  testGraphs(g,gLoad);

  agi::destroyGraph(g);
  agi::destroyGraph(gLoad);
  
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  EnGPar_Finalize();
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}

void testSizes(agi::Ngraph* g1,agi::Ngraph* g2) {
  if (!PCU_Comm_Self())
    printf("Checking Sizes\n");
  assert(g1->numLocalVtxs()==g2->numLocalVtxs());
  assert(g1->numGlobalVtxs()==g2->numGlobalVtxs());
  assert(g1->numLocalEdges(0)==g2->numLocalEdges(0));
  assert(g1->numGlobalEdges(0)==g2->numGlobalEdges(0));
}

void testVertices(agi::Ngraph* g1,agi::Ngraph* g2) {
  if (!PCU_Comm_Self())
    printf("Iterating over vertices\n");
  //Test iterating through vertices
  agi::VertexIterator* gitr1 = g1->begin();
  agi::GraphVertex* vtx1=NULL;
  agi::VertexIterator* gitr2 = g2->begin();
  agi::GraphVertex* vtx2=NULL;
  size_t i=0;
  while ((vtx1 = g1->iterate(gitr1))&&(vtx2 = g2->iterate(gitr2))) {
    i++;
    assert(g1->weight(vtx1)==g2->weight(vtx2));
    assert(g1->degree(vtx1,0)==g2->degree(vtx2,0));
    assert(i<=g1->numLocalVtxs());
  }
  assert(i==g1->numLocalVtxs());

}
void testEdges(agi::Ngraph* g1,agi::Ngraph* g2) {
  if (!PCU_Comm_Self())
    printf("Iterating over edges & pins\n");


  
  //Test iterating through edges & pins on vertices
  agi::EdgeIterator* itr1 = g1->begin(0);
  agi::GraphEdge* edge1=NULL;
  agi::EdgeIterator* itr2 = g2->begin(0);
  agi::GraphEdge* edge2=NULL;
  while ((edge1 = g1->iterate(itr1))&&(edge2 = g2->iterate(itr2))) {
    assert(g1->globalID(edge1)==g2->globalID(edge2));
    assert(g1->weight(edge1)==g2->weight(edge2));
    int np1 = g1->degree(edge1);
    int np2 = g2->degree(edge2);
    assert(np1==np2);
    agi::GraphVertex* out1;
    agi::PinIterator* pitr1 = g1->pins(edge1);
    agi::GraphVertex* out2;
    agi::PinIterator* pitr2 = g2->pins(edge2);
    for (int k=0;k<np1;k++) {
      out1 = g1->iterate(pitr1);
      out2 = g2->iterate(pitr2);
      assert(g1->owner(out1)==g2->owner(out2));
      assert(g1->globalID(out1)==g2->globalID(out2));
    }
  }
}
