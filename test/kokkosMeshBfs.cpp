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
  apf::Field* f = apf::createLagrangeField(m, "vtx_distance", apf::SCALAR, 1);
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
  printf("numEnts %d size of lists %d\n", numEnts, pg->eve_offsets[t][numEnts]);
  LIDs lists_d("lists_d", pg->eve_offsets[t][numEnts]);
  hostToDevice(lists_d, pg->eve_lists[t]);
  //create a device array for the visited flag and initialize it
  LIDs visited_d("visited_d", numEnts);
  LIDs visited_next_d("visited_next_d", numEnts);
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    visited_d(e) = 0;
    if(e==0)
      visited_next_d(e) = 1; //the starting vertex
    else
      visited_next_d(e) = 0;
  });
  //create a device array for the distance and initialize it
  LIDs dist_d("dist_d", numEnts);
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    dist_d(e) = 0;
  });
  //create a flag to indicate if any distances have changed
  agi::lid_t found = 0;
  LIDs found_d("found_d", 1);
  //count iterations
  int iter = 0;
  int max_iter = 1000;
  LIDs iter_d("iter_d", 1);
  //run bfs
  do {
    found = 0; // reset the found flag
    iter++;  // increment the iteration count
    hostToDevice(iter_d, &iter);
    hostToDevice(found_d, &found);
    //
    // insert bfs code here
    Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int u) {
      visited_d(u) |= visited_next_d(u);
      visited_next_d(u) = 0;
    });

    Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int u) {
      // loop over adjacent vertices
      for (int adj = offsets_d(u); adj < offsets_d(u+1); adj++) {
        int v = lists_d(adj);
        visited_next_d(u) |= visited_d(v);
      }
      dist_d(u) = dist_d(u) + ((visited_next_d(u) && !visited_d(u)) * iter_d(0));
      // set the flag to indicate that a new vertex was found 
      if (visited_next_d(u) && !visited_d(u)) found_d(0) = 1;
    });
    deviceToHost(found_d, &found);
    //
  } while (found && iter < max_iter);
  //transfer the distance to the host
  fprintf(stderr, "\n%d total iterations\n", iter);
  agi::lid_t* dist = new int[numEnts];
  deviceToHost(dist_d, dist);
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
  while ((ent = m->iterate(vitr))) {
    int d = distance[i];
    m->setIntTag(ent,dist,&d);
    i++;
  }
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
