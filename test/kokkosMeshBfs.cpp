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
#include <stdio.h>

using engpar::hostToDevice;
using engpar::deviceToHost;
using engpar::LIDs;


// bfs for getting the all depths from starts
void bfs_depth(const LIDs &offsets_d, const LIDs &lists_d, LIDs &dist_d, const LIDs &seeds_d,
    const LIDs &componentIds_d, const LIDs &componentIdStartDepths_d, const agi::lid_t numEnts) {

  const int rank = PCU_Comm_Self();
  
  //create a device array for the visited flag
  LIDs visited_d("visited_d", numEnts);
  LIDs visited_next_d("visited_next_d", numEnts);
  // initialize starts
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) { 
    visited_next_d(e) = seeds_d(e);
  });

  //create a flag to indicate if any distances have changed
  agi::lid_t found = 0;
  LIDs found_d("found_d", 1);
  //count iterations
  int iter = 0;
  const int max_iter = 1000;
  //run bfs
  do {
    found = 0; // reset the found flag
    iter++;  // increment the iteration count
    hostToDevice(found_d, &found);
    //bfs main loop
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
        dist_d(u) = componentIdStartDepths_d(componentIds_d(u)) + iter;
      }
    });
    deviceToHost(found_d, &found);
  } while (found && iter < max_iter);
  if (iter == max_iter) {
    fprintf(stderr, "bfs_depth max iter exceeded\n");
  }
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
  if (iter == max_iter) {
    fprintf(stderr, "bfs_component max iter exceeded\n");
  }
  
  //get size of component
  int size = 0;
  Kokkos::parallel_reduce("SumReduce", numEnts, KOKKOS_LAMBDA(const int i, int& result) {
    result += visited_d(i);
  }, size);
  return size;
}


// functor for getting max outside-in depths of each component
// used with Kokkos::parallel_reduce
class MaxArrayReduceFunctor {
  public:
    typedef int value_type[];
    const unsigned value_count;  // number of components

    const int nwork;  // number of edges
    
    const LIDs componentIds_d;
    const LIDs oiDepths_d;
   
    // nwork: number of edges, count: number of components
    MaxArrayReduceFunctor(const int arg_nwork, const int arg_count, const LIDs &componentIds, const LIDs &oiDepths) : value_count(arg_count), nwork(arg_nwork), componentIds_d(componentIds), oiDepths_d(oiDepths) {}

    KOKKOS_INLINE_FUNCTION
    void init(int dst[]) const {
      for (int i = 0; i < value_count; ++i) dst[i] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    void join(volatile int dst[], const volatile int src[]) const {
      for (int i = 0; i < value_count; ++i) {
        dst[i] = src[i] > dst[i] ? src[i] : dst[i];
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int e, int dst[]) const {
      const int compId = componentIds_d(e);
      const int depth = oiDepths_d(e);
      dst[compId] = depth > dst[compId] ? depth : dst[compId];
    }
};


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
  LIDs componentIds_d("componentIds_d", numEnts); // initialize all component ids as -1
  Kokkos::parallel_for(numEnts, KOKKOS_LAMBDA(const int e) {
    componentIds_d(e) = -1;
  });

  agi::lid_t* componentIds = new int[numEnts];
  std::vector<int> componentSizes;
  agi::lid_t startEdge = 0;
  int compId = 0;
  bool unlabeledEdgeExists = true;
  while (unlabeledEdgeExists) { // while any edge is unlabeled
    int component_size = bfs_component(offsets_d, lists_d, componentIds_d, startEdge, numEnts, compId);
    componentSizes.push_back(component_size);
    deviceToHost(componentIds_d, componentIds);
    // find next start edge
    int i;
    for (i = startEdge+1; i < numEnts; i++) {
      if (componentIds[i] == -1) {
        startEdge = i;
        break;
      }
    }
    if (i >= numEnts) {
      unlabeledEdgeExists = false;
    }
    compId++;
  }

  const int numComponents = compId;

  /* print rank and number of components
  int r;
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  fprintf(stderr, "rank: %d, %d components\n", r, numComponents);
  */

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
  LIDs seeds_d("seeds_d", numEnts);
  Kokkos::parallel_for("compDist_setOiDepthStarts", numStarts, KOKKOS_LAMBDA(const int s) {
    oiDepth_d(starts_d(s)) = 0;
    seeds_d(starts_d(s)) = 1;
  });
  LIDs componentIdStartDepths_d("componentIdStartDepths_d", numComponents);
  bfs_depth(offsets_d, lists_d, oiDepth_d, seeds_d, componentIds_d, componentIdStartDepths_d, numEnts);
/*
  // find max outside-in depth for each component
  LIDs maxOiDepths_d("maxOiDepths_d", numComponents);
  int *maxOiDepths = new int[numComponents];
  for (int componentId = 0; componentId < numComponents; ++componentId) {
    int maxOiDepth = -1; 
    Kokkos::parallel_reduce("maxDepthOfComponent", numEnts, KOKKOS_LAMBDA(const int &e, int& max) {
      int d = (componentIds_d(e) == componentId)*oiDepth_d(e);
      max = (d > max) ? d : max;
    }, Kokkos::Max<int>(maxOiDepth));
    maxOiDepths[componentId] = maxOiDepth;
  }
*/

  // find max outside-in depth for each component
  int *maxOiDepths = new int[numComponents];
  Kokkos::parallel_reduce(numEnts,
                          MaxArrayReduceFunctor(numEnts, numComponents, componentIds_d, oiDepth_d),
                          maxOiDepths);
  LIDs maxOiDepths_d("maxOiDepths_d", numComponents);

  hostToDevice(maxOiDepths_d, maxOiDepths);
  free(maxOiDepths);

  // initialize inside-out depths array and seeds (starts of i-o BFS) array
  // mark edges with max outside-in depth; these are the edges at the
  // components topological center and used as the start of inside-out BFS
  hostToDevice(componentIdStartDepths_d, componentIdStartDepths);
  LIDs ioDepth_d("ioDepth_d", numEnts);
  Kokkos::parallel_for("setSeedsAndDepths", numEnts, KOKKOS_LAMBDA(const int e) {
    const int componentId = componentIds_d(e);
    const bool isMaxDepth = oiDepth_d(e) == maxOiDepths_d(componentId);
    seeds_d(e) = isMaxDepth ? 1 : 0;
    ioDepth_d(e) = isMaxDepth ? componentIdStartDepths_d(componentId) : -1;
  });

  bfs_depth(offsets_d, lists_d, ioDepth_d, seeds_d, componentIds_d, componentIdStartDepths_d, numEnts);

  free(componentIdStartDepths);
  free(componentIds);

  //label each visited hyperedge with its depth in ioDepth
  agi::lid_t* compDist = new agi::lid_t[numEnts];
  deviceToHost(ioDepth_d, compDist);
  return compDist;
}


apf::Field* convert_my_tag(apf::Mesh* m, const char* name, apf::MeshTag* t) {
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = m->begin(0);
  apf::Field* f = apf::createLagrangeField(m, name, apf::SCALAR, 1);
  int vals;
  while ((vtx = m->iterate(it))) {
    m->getIntTag(vtx, t, &vals);
    double v = static_cast<double>(vals);
    apf::setComponents(f, vtx, 0, &v);
  }
  m->end(it);
  return f;
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

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int numRanks;
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  
  double total_time = 0;
  int repeat = 20;
  for (int i = 0; i < repeat; ++i) {
  for (int p = 0; p < numRanks; p++) {
    if (p == rank) {
      Kokkos::Timer timer;
      computeComponentDistance(g, edges[0]);
      double time = timer.seconds();
      total_time += time;
      //fprintf(stderr, "%d %f\n", rank, time);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  }

  for (int p = 0; p < numRanks; ++p) {
    if (p == rank) {
      fprintf(stderr, "%d %f\n", rank, total_time/static_cast<double>(repeat));
    }
  }

  agi::lid_t* compDist = computeComponentDistance(g, edges[0]);
  // tag mesh vertices with component distance
  apf::MeshTag* compDistTag = m->findTag("compDistance");
  if (compDistTag) {
    fprintf(stderr, "compDistance tag found, checking distances\n");
  } else {
    fprintf(stderr, "no compDistance tag found, not checking distances\n");
  }

  apf::MeshTag* engparDistTag = m->createIntTag("engparDistance",1);
  apf::MeshIterator* vitr = m->begin(0);
  apf::MeshEntity* ent;
  int i=0;
  while ((ent = m->iterate(vitr))) {
    const int engDist = compDist[i];
    m->setIntTag(ent,engparDistTag,&engDist);
    if (compDistTag) {
      int d;
      m->getIntTag(ent,compDistTag,&d);
      if(d != compDist[i]) {
        fprintf(stderr,
            "%d component distance does not match reference: %d d%d d_ref%d\n",
            PCU_Comm_Self(), i, compDist[i], d);
      }
    }
    i++;
  }
  m->end(vitr);
  free(compDist);
  convert_my_tag(m,"engparDist",engparDistTag);
  apf::writeVtkFiles("ocean",m);
  destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}
