#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <engpar_kokkosColoring.h>
#include <set>


bool checkDirected(agi::Ngraph* g, agi::etype t=0) {
  agi::PNgraph* pg = g->publicize();
  const agi::lid_t numverts = pg->num_local_verts; 
  const agi::lid_t* degree_list = pg->degree_list[t];
  const agi::lid_t* edge_list = pg->edge_list[t];
  for (int i=0; i<numverts; ++i) { 
    for (int j=degree_list[i]; j<degree_list[i+1]; ++j) {
      bool back_edge = false; 
      for (int k=degree_list[edge_list[j]]; k<degree_list[1+edge_list[j]]; ++k) {
        if (i == edge_list[k])
          back_edge = true;
      }
      if (!back_edge) {
        printf ("Directed graphs are not supported\n");
        return false;
      }
    }
  }
  return true;
}


void vertColor(agi::Ngraph* g, agi::etype t=0) {
  // Call EnGPar graph coloring on vertices
  engpar::ColoringInput* in = engpar::createColoringInput(g, t, true);
  agi::lid_t** colors = new agi::lid_t*[1];
  engpar::EnGPar_KokkosColoring(in, colors); 
  // Assign colors to graph vertices
  agi::GraphTag* tag = g->createIntTag(-1);
  agi::GraphVertex* v;
  agi::VertexIterator* vitr = g->begin();
  int i=0;
  while ((v=g->iterate(vitr)))
    g->setIntTag(tag, v, (*colors)[i++]);

  // Check that the coloring is valid
  agi::GraphEdge* e;
  agi::EdgeIterator* eitr = g->begin(0);
  int conflicts = 0;
  while ((e=g->iterate(eitr))) {
    int u = g->getIntTag(tag, g->u(e));
    int v = g->getIntTag(tag, g->v(e));
    if (u==v) ++conflicts;
  }
  g->destroy(eitr);
  assert(conflicts == 0);
  printf("Vertex Coloring Complete\n");
}


void edgeColor(agi::Ngraph* g, agi::etype t=0) {
  // Call EnGPar graph coloring on edges
  engpar::ColoringInput* in = engpar::createColoringInput(g, t, false);
  agi::lid_t** colors = new agi::lid_t*[1];
  engpar::EnGPar_KokkosColoring(in, colors);
  // assign colors to edges
  agi::GraphTag* tag = g->createIntTag(t);
  agi::GraphEdge* e;
  agi::EdgeIterator* eitr = g->begin(t);
  int i=0;
  while ((e=g->iterate(eitr)))
    g->setIntTag(tag, e, (*colors)[i++]);
  g->destroy(eitr);
  // Check that the coloring is valid
  agi::GraphVertex* v;
  agi::VertexIterator* vitr = g->begin();
  size_t conflicts = 0;
  while ((v=g->iterate(vitr))) {
    agi::EdgeIterator* eitr = g->edges(v, t);
    std::set<int> l_colors;
    for (agi::lid_t i=0; i<g->degree(v, t); ++i) {
      e = g->iterate(eitr);
      l_colors.insert(g->getIntTag(tag, e)); 
    }
    g->destroy(eitr);
    if (l_colors.size() < (size_t) g->degree(v, t))
      ++conflicts;
  }
  assert(conflicts == 0);
  printf("Edge Coloring Complete\n");
}


int main(int argc, char* argv[]) {

  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  if ( argc != 2 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file>",argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  Kokkos::initialize(argc,argv);
  
  agi::Ngraph* g = agi::createBinGraph(argv[1]);
  //assert(checkDirected(g));
  vertColor(g);
  for (agi::lid_t t=0; t<g->numEdgeTypes(); ++t) {
    vertColor(g, t);
    edgeColor(g, t);
  }
  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
