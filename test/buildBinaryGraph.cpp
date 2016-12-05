#include <binGraph.h>
#include <mpi.h>
#include <PCU.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <set>
void testSizes(agi::binGraph* g);
void testVertices(agi::binGraph* g);
void testEdges(agi::binGraph* g);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2&&argc!=3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file> [vertex_partition_file]",argv[0]);
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }
  agi::binGraph* g;
  if (argc==2)
    g = new agi::binGraph(argv[1]);
  else
    g = new agi::binGraph(argv[1],argv[2]);
  testSizes(g);
  testVertices(g);
  testEdges(g);
  delete g;
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}

void testSizes(agi::binGraph* g) {
  if (!PCU_Comm_Self())
    printf("Checking Sizes\n");

  agi::gid_t temp_size = PCU_Add_Long(g->numLocalVtxs());
  assert(temp_size==g->numGlobalVtxs());
  assert(g->numLocalVtxs()+g->numGhostVtxs()<=g->numGlobalVtxs());
  temp_size = PCU_Add_Long(g->numLocalEdges());
  assert(temp_size==g->numGlobalEdges());

}
void testVertices(agi::binGraph* g) {
  if (!PCU_Comm_Self())
    printf("Iterating over vertices\n");
  agi::VertexIterator* gitr = g->begin();
  agi::GraphVertex* vtx=NULL;
  size_t i=0;
  while (vtx = g->iterate(gitr)) {
    i++;
    //assert(g->weight(vtx)==1.0);
    assert(g->degree(vtx,0)>=0);
    assert(i<=g->numLocalVtxs());
  }
  assert(i==g->numLocalVtxs());

}

void testEdges(agi::binGraph* g) {
  agi::VertexIterator* gitr = g->begin();
  agi::GraphVertex* vtx=NULL;
  int num_edges =0;
  int num_ghosts = 0;
  std::set<agi::GraphVertex*> ghosts;
  while (vtx = g->iterate(gitr)) {
    int deg = g->degree(vtx,0);
    agi::GraphVertex* other;
    agi::EdgeIterator* eitr = g->edges(vtx,0);
    for (int j=0;j<deg;j++) {
      other = g->v(g->iterate(eitr));
      if (g->owner(other)!=PCU_Comm_Self())
        ghosts.insert(other);
      
      //assert(g->weight(edge)>0);
      num_edges++;
    }
  }
  assert(ghosts.size()==g->numGhostVtxs());
  assert(num_edges==g->numLocalEdges());
}
