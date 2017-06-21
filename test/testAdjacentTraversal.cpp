#include <apfGraph.h>
#include <apfMesh2.h>
#include <cassert>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <binGraph.h>
#include <apf.h>
#include <vector>
#include <engpar_support.h>

void testAdjacent(agi::Ngraph* g,agi::etype t=0);
void compareTraversal(agi::Ngraph* g,agi::etype t=0);
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if ( argc != 4) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <graph>\n", argv[0]);
    EnGPar_Finalize();
    assert(false);
  }

  //Load in PUMI mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  //Construct Ngraph with edges over mesh faces
  agi::Ngraph* g = agi::createAPFGraph(m,3,2);

  //run the adjacency test
  testAdjacent(g);
  compareTraversal(g);
  
  //Destroy Ngraph
  agi::destroyGraph(g);

  //Construct Ngraph with edges over mesh vertices
  int secondaries[2] = {0,1};
  g = agi::createAPFGraph(m,3,secondaries,2);

  //run the adjacency test
  for (int i=0;i<2;i++) {
    testAdjacent(g,secondaries[i]);
    compareTraversal(g,secondaries[i]);
  }
  
  //Destroy Ngraph
  agi::destroyGraph(g);

  //Destroy Mesh
  m->destroyNative();
  apf::destroyMesh(m);

  //Construct a ngraph from traditional graph
  g = agi::createBinGraph(argv[3]);

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

bool phi(char* word) {
  if (!PCU_Comm_Self())
    printf("hi %s\n", word);
  return true;
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
