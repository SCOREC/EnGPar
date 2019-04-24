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
#include <algorithm>

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

agi::lid_t* bfs(agi::Ngraph* g, agi::lid_t t, std::vector<int> start_vs) {
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
  //create a device array for the distance and initialize it
  LIDs dist_d("dist_d", numEnts);
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    visited_d(e) = 0;
    visited_next_d(e) = 0;
    dist_d(e) = 0;
  });
  //label starting vertices
  agi::lid_t numStarts = static_cast<int>(start_vs.size());
  //convert std::vector to agi::lid_t array for hostToDevice
  agi::lid_t *starts = (agi::lid_t*)calloc(numStarts, sizeof(agi::lid_t));
  for (int i = 0; i < numStarts; i++) {
    starts[i] = start_vs[i];
  }
  LIDs starts_d("starts_d", numStarts);
  hostToDevice(starts_d, starts);
  Kokkos::parallel_for(numStarts, KOKKOS_LAMBDA(const int e) {
    visited_next_d(starts_d(e)) = 1;
  });
  free(starts);
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


// bfs for getting the depths from starts
void bfs_depth(LIDs &offsets_d, LIDs &lists_d, LIDs &dist_d, std::vector<int> start_vs, agi::lid_t numEnts) {
  //create a device array for the visited flag and initialize it
  LIDs visited_d("visited_d", numEnts);
  LIDs visited_next_d("visited_next_d", numEnts);
  //create a device array for the distance and initialize it
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    visited_d(e) = 0;
    visited_next_d(e) = 0;
  });
  //label starting vertices
  agi::lid_t numStarts = static_cast<int>(start_vs.size());
  //convert std::vector to agi::lid_t array for hostToDevice
  agi::lid_t *starts = (agi::lid_t*)calloc(numStarts, sizeof(agi::lid_t));
  for (int i = 0; i < numStarts; i++) {
    starts[i] = start_vs[i];
  }
  LIDs starts_d("starts_d", numStarts);
  hostToDevice(starts_d, starts);
  Kokkos::parallel_for(numStarts, KOKKOS_LAMBDA(const int e) {
    visited_next_d(starts_d(e)) = 1;
  });
  free(starts);
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
}


// bfs for labeling components with a single id
int bfs_component(LIDs &offsets_d, LIDs &lists_d, LIDs &component_ids_d, agi::lid_t start, agi::lid_t numEnts, int component_id) {
  LIDs visited_d("visited_d", numEnts);
  LIDs visited_next_d("visited_next_d", numEnts);
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    visited_d(e) = 0;
    if (e == start) {
      visited_next_d(e) = 1;
    } else {
      visited_next_d(e) = 0;
    }
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
    // bfs main loop
    Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int u) {
      visited_d(u) |= visited_next_d(u);
      visited_next_d(u) = 0;
      // loop over adjacent vertices
      for (int adj = offsets_d(u); adj < offsets_d(u+1); adj++) {
        int v = lists_d(adj);
        visited_next_d(u) |= visited_d(v);
      }
      component_ids_d(u) = component_id;
      // set the flag to indicate that a new vertex was found 
      if (visited_next_d(u) && !visited_d(u)) found_d(0) = 1;
    });
    deviceToHost(found_d, &found);
  } while (found && iter < max_iter);
  fprintf(stderr, "\n%d total iterations\n", iter);
  //get size of component
  int size = 0;
  Kokkos::parallel_reduce( "SumReduce", numEnts, KOKKOS_LAMBDA(const int i, int& result) {
    result += visited_d(i);
  }, size);
  fprintf(stderr, "\ncomponent size: %d\n", size);
  return size;

}



void computeComponentDistance(agi::Ngraph* g, agi::lid_t t) {
  // graph setup
  assert(g->isHyper());
  agi::PNgraph* pg = g->publicize();
  g->parallel_create_eve(t);
  agi::lid_t numEnts = pg->num_local_edges[t]; // size
  //transfer the eve to the device
  LIDs offsets_d("offsets_d", numEnts+1);
  hostToDevice(offsets_d, pg->eve_offsets[t]);
  printf("numEnts %d size of lists %d\n", numEnts, pg->eve_offsets[t][numEnts]);
  LIDs lists_d("lists_d", pg->eve_offsets[t][numEnts]);
  hostToDevice(lists_d, pg->eve_lists[t]);

  // determine component for each hyperedge
  std::vector<int> componentSizes;
  LIDs componentIds_d("componentIds_d", numEnts); // initialize all component ids as -1
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    componentIds_d(e) = -1;
  });
  agi::lid_t* component_ids = new int[numEnts];
  agi::lid_t startEdge = 0;
  agi::lid_t prevStart = 0;
  int compId = 0;
  while (startEdge != -1) { // while any edge is unlabeled
    int component_size = bfs_component(offsets_d, lists_d, componentIds_d, startEdge, numEnts, compId);
    componentSizes.push_back(component_size);
    deviceToHost(componentIds_d, component_ids);
    // find next start edge
    startEdge = -1;
    for (int i = prevStart; i < numEnts; i++) {
      if (component_ids[i] == -1) {
        startEdge = i;
        break;
      }
    }
    prevStart = startEdge; 
    compId++;
  }

  int numComponents = compId;
  printf("%d components found\n", numComponents);
  for (int i = 0; i < numComponents; i++) {
    printf("component %d: size %d\n", i, componentSizes[i]);
  }

  // compute starting depth for each component on the host
  // sort componentSizes array
  std::vector<int> componentIndices; 
  for (int i = 0; i < numComponents; i++) {
    componentIndices.push_back(i);
  }
  std::sort(componentIndices.begin(), componentIndices.end(), [&](int i, int j){return componentSizes[i]>componentSizes[j];});
  std::sort(componentSizes.begin(), componentSizes.end(), std::greater<int>());
  // get starting depths
  int *startingDepth = new int[numComponents];  // descending exclusive sum values
  int *componentIdStartDepths = new int[numComponents];  // index i is startingDepth of component i
  int sizeSum = componentSizes[0];
  for (int i = 0; i < numComponents; i++) {
    printf("%d\n", sizeSum);
    startingDepth[i] = sizeSum;
    sizeSum += componentSizes[i];
    componentIdStartDepths[componentIndices[i]] = startingDepth[i];
  } 
  free(startingDepth); 
 
  
  // run the outside-in BFS to find the core vertices
  std::vector<int> startEdges;
  for (int i = 0; i < numEnts; i++) {
    if (g->isCut(pg->getEdge(i, t))) startEdges.push_back(i);
  }
  printf("\n%u start edges for outside in\n", startEdges.size());
  LIDs oiDepth_d("oiDepth_d", numEnts);
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    oiDepth_d(e) = -1;
  });

  // set the depth of each edge in startEdges to 0 in oiDepth
  agi::lid_t numStarts = static_cast<int>(startEdges.size());
  //convert std::vector to agi::lid_t array for hostToDevice
  agi::lid_t *starts = (agi::lid_t*)calloc(numStarts, sizeof(agi::lid_t));
  for (int i = 0; i < numStarts; i++) {
    starts[i] = startEdges[i];
  }
  LIDs starts_d("starts_d", numStarts);
  hostToDevice(starts_d, starts);
  Kokkos::parallel_for(numStarts, KOKKOS_LAMBDA(const int e) {
    oiDepth_d(starts_d(e)) = 0;
  });
  free(starts);
  bfs_depth(offsets_d, lists_d, oiDepth_d, startEdges, numEnts);
  
  // run the inside-out BFS to compute the distance from the component center to the boundary
  LIDs ioDepth_d("ioDepth_d", numEnts);
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    ioDepth_d(e) = -1;
  });
  LIDs componentIdStartDepths_d("componentIdStartDepths_d", numComponents);
  hostToDevice(componentIdStartDepths_d, componentIdStartDepths);
  LIDs compId_d("compId_d", 1);
  LIDs maxVal_d("maxVal_d", 1);
  for (int compId = 0; compId < numComponents; compId++) {
    // host compId to device
    hostToDevice(compId_d, &compId);
    // parallel_reduce
    int maxVal;
    Kokkos::parallel_reduce("MaxReduce", numEnts, KOKKOS_LAMBDA(const int& e, int& max) {
      int d = (componentIds_d(e) == compId_d(0))*oiDepth_d(e);
      if (d > max) max = d;
    }, Kokkos::Max<int>(maxVal));
    // host maxVal to device
    hostToDevice(maxVal_d, &maxVal);
    // parallel_for
    LIDs seeds_d("seeds_d", numEnts);
    Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
      seeds_d(e) = (oiDepth_d(e) == maxVal_d(0))*(componentIds_d(e) == compId_d(0));
      ioDepth_d(e) = (oiDepth_d(e) == maxVal_d(0))*(componentIds_d(e) == compId_d(0))*componentIdStartDepths_d(compId_d(0));
    });
    // get seeds for io BFS
    int *seeds = new int[numEnts];
    deviceToHost(seeds_d, seeds);
    std::vector<int> startEdges;
    for (int i = 0; i < numEnts; i++) {
      if (seeds[i]) startEdges.push_back(i);
    }
    // run bfs
    bfs_depth(offsets_d, lists_d, ioDepth_d, startEdges, numEnts);
  }
  // bfs label each visited hyperedge with its depth in ioDepth
  
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
  computeComponentDistance(g, edges[0]);
  // run bfs
  agi::lid_t* distance = bfs(g,edges[0],{ 0, 20, 40 });
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
