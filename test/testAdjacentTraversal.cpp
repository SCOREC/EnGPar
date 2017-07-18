#include <cassert>
#include <binGraph.h>
#include <vector>
#include <engpar_support.h>

void testAdjacent(agi::Ngraph* g,agi::etype t=0);
void compareTraversal(agi::Ngraph* g,agi::etype t=0);
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if ( argc != 3) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <bgd_prefix> <graph>\n", argv[0]);
    EnGPar_Finalize();
    assert(false);
  }

  //Construct Ngraph with edges over mesh faces
  agi::Ngraph* g = agi::createEmptyGraph();
  g->loadFromFile(argv[1]);

  //run the adjacency test
  testAdjacent(g);
  compareTraversal(g);
  
  //Destroy Ngraph
  agi::destroyGraph(g);

  //Construct a ngraph from traditional graph
  g = agi::createBinGraph(argv[2]);

  //run the adjacency test
  testAdjacent(g);
  compareTraversal(g);

  //Destroy the graph
  agi::destroyGraph(g);
  
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}

void testAdjacent(agi::Ngraph* g,agi::etype t) {
  if (!PCU_Comm_Self())
    printf("Beginning Traversal\n");
  agi::VertexIterator* vitr = g->begin();
  agi::GraphVertex* vtx = NULL;
  agi::lid_t num_pins =0;
  agi::lid_t num_edges=0;
  while ((vtx = g->iterate(vitr))) {
    agi::GraphIterator* gitr = g->adjacent(vtx,t);
    agi::GraphVertex* other = NULL;
    agi::GraphEdge* edge = NULL;
    while ((other = g->iterate(gitr))) {
      assert(other);
      edge = g->edge(gitr);
      assert(edge);
      if (g->isEqual(vtx,other))
        num_pins++;
      num_edges++;      
    }
    g->destroy(gitr);
  }
  if (g->isHyper()) {
    if (PCU_Comm_Peers()==1)
      assert(num_pins==g->numLocalPins(t));
  }
  else {
    assert(num_edges==g->numLocalEdges(t));
  }
}

void compareTraversal(agi::Ngraph* g,agi::etype t) {
  if (!PCU_Comm_Self())
    printf("Beginning Comparison\n");
  agi::VertexIterator* vitr = g->begin();
  agi::GraphVertex* vtx = NULL;
  std::vector<agi::GraphEdge*> edges;
  std::vector<agi::GraphVertex*> vtxs;
  while ((vtx = g->iterate(vitr))) {
    agi::GraphIterator* gitr = g->adjacent(vtx,t);
    agi::GraphVertex* other = NULL;
    while ((other = g->iterate(gitr))) {
      edges.push_back(g->edge(gitr));
      vtxs.push_back(other);
    }
    g->destroy(gitr);
  }
  vitr = g->begin();
  agi::lid_t i=0;
  while ((vtx = g->iterate(vitr))) {
    agi::EdgeIterator* eitr = g->edges(vtx,t);
    agi::GraphVertex* other;
    agi::GraphEdge* edge;
    while ((edge = g->iterate(eitr))) {
      assert(edges[i]==edge);
      if (g->isHyper()) {
        agi::PinIterator* pitr = g->pins(edge);
        for (agi::lid_t j=0;j<g->degree(edge);j++) {
          other = g->iterate(pitr);
          assert(g->isEqual(other,vtxs[i]));
          assert(edges[i]==edge);
          i++;
        }
        g->destroy(pitr);
      }
      else {
        assert(g->v(edge)==vtxs[i]);
        i++;
      }
    }
    g->destroy(eitr);
  }
  assert(i==vtxs.size());
  assert(i==edges.size());
}
