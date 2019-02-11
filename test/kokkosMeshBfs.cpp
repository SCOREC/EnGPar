#include <engpar_support.h>
#include <apfGraph.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <cstring>
#include <set>

using engpar::hostToDevice;
using engpar::deviceToHost;
using engpar::LIDs;

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

agi::lid_t* bfs(agi::Ngraph* g, agi::lid_t t) {
  assert(g->isHyper());
  agi::PNgraph* pg = g->publicize();
  g->parallel_create_eve(t);
  agi::lid_t numEnts = pg->num_local_edges[t];
  //transfer the eve to the device
  LIDs offsets_d("offsets_d", numEnts+1);
  hostToDevice(offsets_d, pg->eve_offsets[t]);
  LIDs lists_d("lists_d", pg->eve_lists[t][numEnts]);
  hostToDevice(lists_d, pg->eve_lists[t]);
  //create a device array for the distance
  LIDs dist_view("dist_view", numEnts);
  //write something into the distance - replace this with pull bfs
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
      dist_view(e) = e;
  });
  //transfer the distance to the host
  agi::lid_t* dist = new agi::lid_t[numEnts];
  deviceToHost(dist_view, dist);
  return dist;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  Kokkos::initialize(argc,argv);
  if ( argc != 3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh>",argv[0]);
    Kokkos::finalize();
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }
  // Load mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]); 
  int edges[1] = {0}; 
  agi::Ngraph* g = agi::createAPFGraph(m,"mesh_graph",3,edges,1);
  // run bfs
  agi::lid_t* distance = bfs(g,edges[0]);
  // Move colors to mesh
  apf::MeshTag* dist = m->createIntTag("distance",1);
  apf::MeshIterator* vitr = m->begin(0);
  apf::MeshEntity* ent;
  int i=0;
  while ((ent = m->iterate(vitr)))
    m->setIntTag(ent,dist,&distance[i++]);
  m->end(vitr);
  convert_my_tag(m, dist);
  apf::writeVtkFiles("meshDistance", m);

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}
