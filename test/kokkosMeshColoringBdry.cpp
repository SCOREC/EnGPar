#include <engpar_support.h>
#include <apfGraph.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <pcu_util.h> // provides PCU_ALWAYS_ASSERT
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

void checkColoring(apf::Mesh* m, apf::MeshTag* c) {
  apf::MeshIterator* vtxItr = m->begin(0);
  apf::MeshEntity* u;
  while ((u = m->iterate(vtxItr))) {
    if( ! m->isShared(u) )
      continue;
    if( ! m->isOwned(u) )
      continue;
    int uColor;
    m->getIntTag(u,c,&uColor);
    apf::Adjacent adjVerts;
    apf::getBridgeAdjacent(m,u,m->getDimension(),0,adjVerts);
    for(size_t i=0; i<adjVerts.size(); i++) {
      apf::MeshEntity* v = adjVerts[i];
      if( m->isShared(v) && m->isOwned(v) ) {
        int vColor;
        m->getIntTag(v,c,&vColor);
        PCU_ALWAYS_ASSERT(uColor != vColor);
      }
    }
  }
  m->end(vtxItr);
}


agi::lid_t* edgeColor(agi::Ngraph* g, agi::etype t=0) {
  // Call EnGPar graph coloring on edges
  engpar::ColoringInput* in = engpar::createBdryColoringInput(g, t);
  agi::lid_t** colors = new agi::lid_t*[1];
  engpar::EnGPar_KokkosColoring(in, colors);
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
  agi::Ngraph* g = agi::createAPFGraph(m,"mesh_graph",m->getDimension(),edges,1);
  agi::lid_t** colors = new agi::lid_t*[g->numEdgeTypes()];
  for (agi::lid_t t=0; t<g->numEdgeTypes(); ++t) {
    colors[t] = edgeColor(g, t); 
  }
  // Move colors to mesh
  apf::MeshTag* coloring = m->createIntTag("coloring",1);
  apf::MeshIterator* vitr = m->begin(0);
  apf::MeshEntity* ent;
  int i=0;
  while ((ent = m->iterate(vitr))) {
    m->setIntTag(ent,coloring,&colors[0][i++]);
  }
  m->end(vitr);
  checkColoring(m,coloring);
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
