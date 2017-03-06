#include <binGraph.h>
#include <mpi.h>
#include <PCU.h>
#include <cassert>

void traverseEdges(agi::Ngraph*);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2&&argc!=3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file>",argv[0]);
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }

  agi::Ngraph* g = agi::createBinGraph(argv[1]);

  traverseEdges(g);
  
  destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("All tests passed\n");

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}

void traverseEdges(agi::Ngraph* g) {
  agi::EdgeIterator* eitr = g->begin(0);
  agi::GraphEdge* edge =NULL;
  int numEdges=0;
  std::map<agi::GraphVertex*,int> edgesPerVertex;
  while ((edge = g->iterate(eitr))) {
    numEdges++;
    agi::GraphVertex* u = g->u(edge);
    edgesPerVertex[u]++;
    agi::GraphVertex* v = g->v(edge);
    agi::EdgeIterator* eitr2 = g->edges(u);
    agi::GraphEdge* e2 = NULL;
    while ((e2 = g->iterate(eitr2))) {
      if (e2==edge) {
	assert(g->isEqual(v,g->v(e2)));
	assert(g->isEqual(u,g->u(e2)));
      }
    }
  }
  assert(numEdges==g->numLocalEdges());
  std::map<agi::GraphVertex*,int>::iterator itr;
  for (itr=edgesPerVertex.begin();itr!=edgesPerVertex.end();itr++) {
    assert(itr->second==g->degree(itr->first));
  }
  g->destroy(eitr);
}
