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

// bfs for getting the depths from starts
void bfs_depth(LIDs &offsets_d, LIDs &lists_d, LIDs &dist_d, std::vector<int> start_vs,
    agi::lid_t numEnts, agi::lid_t startDepth) {
  const int rank = PCU_Comm_Self();
  //create a device array for the visited flag
  LIDs visited_d("visited_d", numEnts);
  LIDs visited_next_d("visited_next_d", numEnts);
  //label starting vertices
  agi::lid_t numStarts = static_cast<int>(start_vs.size());
  LIDs starts_d("starts_d", numStarts);
  hostToDevice(starts_d, &start_vs[0]);
  Kokkos::parallel_for(numStarts, KOKKOS_LAMBDA(const int e) {
    visited_next_d(starts_d(e)) = 1;
  });
  //create a flag to indicate if any distances have changed
  agi::lid_t found = 0;
  LIDs found_d("found_d", 1);
  //count iterations
  int iter = 0;
  int max_iter = 1000;
  //run bfs
  do {
    found = 0; // reset the found flag
    iter++;  // increment the iteration count
    hostToDevice(found_d, &found);
    // bfs
    Kokkos::parallel_for("bfsDepth_setNext", numEnts, KOKKOS_LAMBDA(const int u) {
      visited_d(u) |= visited_next_d(u);
      visited_next_d(u) = 0;
    });
    Kokkos::parallel_for("bfsDepth_setDist", numEnts, KOKKOS_LAMBDA(const int u) {
      // loop over adjacent vertices
      for (int adj = offsets_d(u); adj < offsets_d(u+1); adj++) {
        int v = lists_d(adj);
        visited_next_d(u) |= visited_d(v);
      }
      // a new vertex
      if( visited_next_d(u) && !visited_d(u) ) {
        found_d(0) = 1;
        dist_d(u) = startDepth+iter;
      }
    });
    deviceToHost(found_d, &found);
  } while (found && iter < max_iter);
}


// bfs for labeling components with a single id
int bfs_component(LIDs &offsets_d, LIDs &lists_d, LIDs &component_ids_d, agi::lid_t start, agi::lid_t numEnts, int component_id) {
  const int rank = PCU_Comm_Self();
  LIDs visited_d("visited_d", numEnts);
  LIDs visited_next_d("visited_next_d", numEnts);
  Kokkos::parallel_for("bfsComp_markStart", numEnts, KOKKOS_LAMBDA(const int e) {
    visited_d(e) = (e == start) ? 1 : 0;
    component_ids_d(e) = (e == start) ? component_id : component_ids_d(e);
  }); 
  //create a flag to indicate if any distances have changed
  agi::lid_t found = 0;
  LIDs found_d("found_d", 1);
  //count iterations
  int iter = 0;
  int max_iter = 1000;
  //run bfs
  do {
    found = 0; // reset the found flag
    iter++;  // increment the iteration count
    hostToDevice(found_d, &found);
    // bfs main loop
    Kokkos::parallel_for("bfsComp_setVisited", numEnts, KOKKOS_LAMBDA(const int u) {
      visited_d(u) |= visited_next_d(u);
      visited_next_d(u) = 0;
    });
    Kokkos::parallel_for("bfsComp_markNext", numEnts, KOKKOS_LAMBDA(const int u) {
      // loop over adjacent vertices
      for (int adj = offsets_d(u); adj < offsets_d(u+1); adj++) {
        int v = lists_d(adj);
        visited_next_d(u) |= visited_d(v);
      }
      // set the flag to indicate that a new vertex was found 
      if (visited_next_d(u) && !visited_d(u)) {
        component_ids_d(u) = component_id;
        found_d(0) = 1;
      }
    });
    deviceToHost(found_d, &found);
  } while (found && iter < max_iter);
  //get size of component
  int size = 0;
  Kokkos::parallel_reduce("SumReduce", numEnts, KOKKOS_LAMBDA(const int i, int& result) {
    result += visited_d(i);
  }, size);
  return size;
}



agi::lid_t* computeComponentDistance(agi::Ngraph* g, agi::lid_t t) {
  const int rank = PCU_Comm_Self();
  // graph setup
  assert(g->isHyper());
  agi::PNgraph* pg = g->publicize();
  g->parallel_create_eve(t);
  agi::lid_t numEnts = pg->num_local_edges[t]; // size
  //transfer the eve to the device
  LIDs offsets_d("offsets_d", numEnts+1);
  hostToDevice(offsets_d, pg->eve_offsets[t]);
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
  int compId = 0;
  while (startEdge != -1) { // while any edge is unlabeled
    int component_size = bfs_component(offsets_d, lists_d, componentIds_d, startEdge, numEnts, compId);
    componentSizes.push_back(component_size);
    deviceToHost(componentIds_d, component_ids);
    // find next start edge
    startEdge = -1;
    for (int i = 0; i < numEnts; i++) {
      if (component_ids[i] == -1) {
        startEdge = i;
        break;
      }
    }
    compId++;
  }

  const int numComponents = compId;
  // compute starting depth for each component on the host
  // sort componentSizes array
  std::vector<int> componentIndices; 
  for (int i = 0; i < numComponents; i++) {
    componentIndices.push_back(i);
  }
  std::sort(componentIndices.begin(), componentIndices.end(), [&](int i, int j){return componentSizes[i]>componentSizes[j];});
  std::sort(componentSizes.begin(), componentSizes.end(), std::greater<int>());
  // get starting depths
  int *componentIdStartDepths = new int[numComponents];  // index i is startingDepth of component i
  int sizeSum = 0;
  for (int i = 0; i < numComponents; i++) {
    componentIdStartDepths[componentIndices[i]] = sizeSum;
    sizeSum += componentSizes[i];
  } 
  
  // run the outside-in BFS to find the core vertices
  std::vector<int> startEdges;
  for (int i = 0; i < numEnts; i++) {
    if (g->isCut(pg->getEdge(i, t))) {
      startEdges.push_back(i);
    }
  }
  LIDs oiDepth_d("oiDepth_d", numEnts);
  Kokkos::parallel_for("compDist_initOiDepth", numEnts, KOKKOS_LAMBDA(const int e) {
    oiDepth_d(e) = -1;
  });

  // set the depth of each edge in startEdges to 0 in oiDepth
  const int numStarts = static_cast<int>(startEdges.size());
  LIDs starts_d("starts_d", numStarts);
  hostToDevice(starts_d, &startEdges[0]);
  Kokkos::parallel_for("compDist_setOiDepthStarts", numStarts, KOKKOS_LAMBDA(const int e) {
    oiDepth_d(starts_d(e)) = 0;
  });
  bfs_depth(offsets_d, lists_d, oiDepth_d, startEdges, numEnts, 0);

  // run the inside-out BFS to compute the distance from the component center to the boundary
  LIDs ioDepth_d("ioDepth_d", numEnts);
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    ioDepth_d(e) = -1;
  });
  LIDs componentIdStartDepths_d("componentIdStartDepths_d", numComponents);
  hostToDevice(componentIdStartDepths_d, componentIdStartDepths);
  for (int compId = 0; compId < numComponents; compId++) {
    // find max outside-in depth for the component
    int maxOiDepth;
    Kokkos::parallel_reduce("maxDepthOfComponent", numEnts, KOKKOS_LAMBDA(const int& e, int& max) {
      int d = (componentIds_d(e) == compId)*oiDepth_d(e);
      max = (d > max) ? d : max;
    }, Kokkos::Max<int>(maxOiDepth));
    // mark the edges with max outside-in depth; these are the edges at the
    // components topological center and are used as the start for the
    // inside-out BFS
    LIDs seeds_d("seeds_d", numEnts);
    Kokkos::parallel_for("setSeedsAndDepth", numEnts, KOKKOS_LAMBDA(const int e) {
      seeds_d(e) = (oiDepth_d(e) == maxOiDepth) && (componentIds_d(e) == compId) ? 1 : 0;
      ioDepth_d(e) = (oiDepth_d(e) == maxOiDepth) && (componentIds_d(e) == compId) ? 
                     componentIdStartDepths_d(compId) : ioDepth_d(e);
    });
    // get seeds for io BFS
    int *seeds = new int[numEnts];
    deviceToHost(seeds_d, seeds);
    std::vector<int> startEdges;
    for (int i = 0; i < numEnts; i++) {
      if (seeds[i]) {
        startEdges.push_back(i);
      }
    }
    // run bfs
    bfs_depth(offsets_d, lists_d, ioDepth_d, startEdges, numEnts, componentIdStartDepths[compId]);
  }
  //label each visited hyperedge with its depth in ioDepth
  agi::lid_t* compDist = new agi::lid_t[numEnts];
  deviceToHost(ioDepth_d, compDist);
  return compDist;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  Kokkos::initialize(argc,argv);
  if ( argc != 3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <in: model> <in: mesh>",argv[0]);
    Kokkos::finalize();
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }
  // Load mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]); 
  const int elmDim = m->getDimension();
  int edges[1] = {0}; 
  agi::Ngraph* g = agi::createAPFGraph(m,"mesh_graph",elmDim,edges,1);
  agi::lid_t* compDist = computeComponentDistance(g, edges[0]);

  // tag mesh vertices with component distance
  apf::MeshTag* compDistTag = m->findTag("compDistance");
  apf::MeshIterator* vitr = m->begin(0);
  apf::MeshEntity* ent;
  int i=0;
  while ((ent = m->iterate(vitr))) {
    int d;
    m->getIntTag(ent,compDistTag,&d);
    if(d != compDist[i]) {
      fprintf(stderr,
          "%d component distance does not match reference: %d d d_ref\n",
          PCU_Comm_Self(), i, compDist[i], d);
      exit(EXIT_FAILURE);
    }
    i++;
  }
  m->end(vitr);
  free(compDist);

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}
