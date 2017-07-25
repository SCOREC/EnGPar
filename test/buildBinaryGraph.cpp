#include <binGraph.h>
#include <mpi.h>
#include <PCU.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
#include <engpar_support.h>
void testSizes(agi::Ngraph* g);
void testVertices(agi::Ngraph* g);
void testEdges(agi::Ngraph* g);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  if ( argc != 2&&argc!=3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file> [vertex_partition_file]",argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  //Construct Graph from binary file
  agi::Ngraph* g;
  if (argc==2)
    g = agi::createBinGraph(argv[1]);
  else
    g = agi::createBinGraph(argv[1],argv[2]);
  
  //Test different members of graph
  testSizes(g);
  testVertices(g);
  testEdges(g);

  //Check validity
  agi::checkValidity(g);
  
  //Destroy the graph
  destroyGraph(g);
  
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}

void testSizes(agi::Ngraph* g) {
  if (!PCU_Comm_Self())
    printf("Checking Sizes\n");

  agi::gid_t temp_size = PCU_Add_Long(g->numLocalVtxs());
  assert(temp_size==g->numGlobalVtxs());
  assert(g->numLocalVtxs()+g->numGhostVtxs()<=g->numGlobalVtxs());
  temp_size = PCU_Add_Long(g->numLocalEdges());
  assert(temp_size==g->numGlobalEdges());

}
void testVertices(agi::Ngraph* g) {
  if (!PCU_Comm_Self())
    printf("Iterating over vertices\n");
  agi::VertexIterator* gitr = g->begin();
  agi::GraphVertex* vtx=NULL;
  int i=0;
  while ((vtx = g->iterate(gitr))) {
    i++;
    //assert(g->weight(vtx)==1.0);
    assert(i<=g->numLocalVtxs());
  }
  assert(i==g->numLocalVtxs());

}

void testEdges(agi::Ngraph* g) {
  if (!PCU_Comm_Self())
    printf("Iterating over edges\n");
  agi::VertexIterator* gitr = g->begin();
  agi::GraphVertex* vtx=NULL;
  agi::lid_t num_edges =0;
  std::set<agi::GraphVertex*> ghosts;
  while ((vtx = g->iterate(gitr))) {
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
    g->destroy(eitr);
    int temp_count =0;
    eitr = g->edges(vtx,0);
    agi::GraphEdge* edge;
    while ((edge = g->iterate(eitr))) {
      other = g->v(edge);
      temp_count++;
    }
    g->destroy(eitr);
    assert(temp_count==deg);
  }
  assert((int)ghosts.size()==g->numGhostVtxs());
  assert(num_edges==g->numLocalEdges());

}

