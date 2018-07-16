#include <engpar_support.h>
#include <apfGraph.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <engpar_kokkosColoring.h>
#include <cstring>
#include <set>


apf::Field* convert_my_tag(apf::Mesh* m, apf::MeshTag* t) {
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = m->begin(0);
  apf::Field* f = apf::createLagrangeField(m, "my_field", apf::SCALAR, 1);
  int val;
  while ((vtx = m->iterate(it))) {
    m->getIntTag(vtx, t, &val);
    apf::setScalar(f, vtx, 0, (double) val);
  }
  m->end(it);
  return f;
}


agi::lid_t* edgeColor(agi::Ngraph* g, agi::etype t=0) {
  // Call EnGPar graph coloring on edges
  engpar::ColoringInput* in = engpar::createColoringInput(g, t);
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
  printf("%lu\n",conflicts);
  return *colors;
}


int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  if ( argc != 3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh>",argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }
  Kokkos::initialize(argc,argv);
  // Load mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]); 
  int edges[1] = {0};
  agi::Ngraph* g = agi::createAPFGraph(m,3,edges,1);
  agi::lid_t** colors = new agi::lid_t*[g->numEdgeTypes()];
  for (agi::lid_t t=0; t<g->numEdgeTypes(); ++t)
    colors[t] = edgeColor(g, t);

  // Move colors to mesh 
  apf::MeshTag* coloring = m->createIntTag("coloring",1);
  apf::MeshIterator* vitr = m->begin(0);
  apf::MeshEntity* ent;
  int i=0;
  while ((ent = m->iterate(vitr))) {
    m->setIntTag(ent,coloring,&colors[0][i++]);
  }
  m->end(vitr);
  convert_my_tag(m, coloring);
  apf::writeVtkFiles("meshColoring", m);


  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}
