#include <engpar_support.h>
#include <engpar.h>
#include <PCU.h>
#include <binGraph.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMesh2.h>
#include <apfGraph.h>
#include <apfMDS.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include <utility>
typedef unsigned int uint;
typedef Kokkos::View<uint*> kkLidView;
typedef Kokkos::TeamPolicy<>::member_type member_type;



void hostToDevice(kkLidView d, agi::lid_t* h) {
  kkLidView::HostMirror hv = Kokkos::create_mirror_view(d);
  for (size_t i=0; i<hv.size(); ++i)
    hv(i) = h[i];
  Kokkos::deep_copy(d,hv);
}


void deviceToHost(kkLidView d, agi::lid_t* h) {
  kkLidView::HostMirror hv = Kokkos::create_mirror_view(d);
  Kokkos::deep_copy(hv,d);
  for(size_t i=0; i<hv.size(); ++i)
    h[i] = hv(i);
}


void parallel_eve(agi::Ngraph* g, agi::etype t=0) {
  agi::PNgraph* pg = g->publicize();
  if (pg->eve_offsets[t]) {
    delete [] pg->eve_offsets[t];
    delete [] pg->eve_lists[t];
  }
  // Load graph info to device
  const int N = pg->num_local_edges[t];
  const int M = pg->num_local_verts; 

  Kokkos::View<uint**> A ("adjacency_matrix", N, N);
  kkLidView degree_view ("degree_view", M+1);
  hostToDevice(degree_view, pg->degree_list[t]); 
  kkLidView edge_view ("edge_view", pg->num_local_pins[t]);
  if (pg->isHyperGraph) {
    hostToDevice(edge_view, pg->edge_list[t]);
    // Uses hierarchical parallelism to construct adjacency matrix 
    // HYPERGRAPH VERSION
    Kokkos::parallel_for (Kokkos::TeamPolicy<>(M, Kokkos::AUTO()),
      KOKKOS_LAMBDA (const member_type& thread) {
        const int v = thread.league_rank();
        const int loop_count = degree_view(v+1) - degree_view(v); 
        Kokkos::parallel_for (Kokkos::TeamThreadRange(thread, loop_count), 
          [=] (const int i) {
            int row = degree_view(v) + i;
            Kokkos::parallel_for (Kokkos::ThreadVectorRange(thread, loop_count),
              [=] (const int j) { 
                int col = degree_view(v) + j;
                if (edge_view(row)!=edge_view(col))
                  A(edge_view(row),edge_view(col)) = 1;
              });
          });
      });
  } else {
    // GRAPH VERSION
    Kokkos::parallel_for (Kokkos::TeamPolicy<>(M, Kokkos::AUTO()),
      KOKKOS_LAMBDA (const member_type& thread) {
        const int v = thread.league_rank();
        const int loop_count = degree_view(v+1) - degree_view(v); 
        Kokkos::parallel_for (Kokkos::TeamThreadRange(thread, loop_count), 
          [=] (const int i) {
            int row = degree_view(v) + i;
            Kokkos::parallel_for (Kokkos::ThreadVectorRange(thread, loop_count),
              [=] (const int j) {
                int col = degree_view(v) + j;
                if (row!=col)
                  A(row,col) = 1;
              });
          });
      });
  }
  // Compress 'A' into a CSR
  kkLidView eve_deg ("A_deg", N+1); 
  Kokkos::parallel_for (Kokkos::TeamPolicy<>(N, Kokkos::AUTO()),
    KOKKOS_LAMBDA (const member_type& thread) {
      const int i = thread.league_rank();
      int degree = 0;
      Kokkos::parallel_reduce (Kokkos::TeamThreadRange(thread, N),
        [=] (const int j, int& ldegree) {
          ldegree += A(i,j);
        }, degree);
      eve_deg(i+1) = degree;
    });
  Kokkos::parallel_scan(N, KOKKOS_LAMBDA(const int& i, int& upd, const bool& final) {
    const int val = eve_deg(i+1); 
    upd += val;
    if (final)
      eve_deg(i+1) = upd;
  });
  pg->eve_offsets[t] = new int[N+1];
  deviceToHost(eve_deg, pg->eve_offsets[t]);
  int eve_edge_size = pg->eve_offsets[t][N]; 
  kkLidView eve_edges ("eve_edges", eve_edge_size);
  Kokkos::parallel_for (N, KOKKOS_LAMBDA(const int i) {
    int index = eve_deg(i);
    while (index != eve_deg(i+1)) {
      for (int j=0; j<N; ++j) {
        if (A(i,j) == 1)
          eve_edges(index++) = j;
      }
    }
  }); 
  pg->eve_lists[t] = new int[eve_edge_size];
  deviceToHost(eve_edges, pg->eve_lists[t]);
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
  // Create graphs 
  agi::Ngraph* g_parallel = agi::createAPFGraph(m,"g_parallel",3,edges,1);
  agi::PNgraph* gp = g_parallel->publicize();
  agi::Ngraph* g_serial = agi::createAPFGraph(m,"g_serial",3,edges,1);
  agi::PNgraph* gs = g_serial->publicize();

  // Modify for multiple edge types 
  // Time the two constructions
  double t0 = PCU_Time();
  parallel_eve(g_parallel);
  printf("Parallel Construction Time: %f\n",PCU_Time()-t0);
  t0 = PCU_Time();
  g_serial->create_eve_adjacency(0);
  printf("Serial Construction Time: %f\n",PCU_Time()-t0);
  // Compare the two constructions 
  int conflicts = 0;
  for (int i=0; i<gp->num_local_edges[0]+1; ++i) {
    if (gp->eve_offsets[0][i] != gs->eve_offsets[0][i])
      ++conflicts;
  }
  printf("Offsets conflicts: %i\n", conflicts);
  conflicts = 0;
  for (int e=0; e<gs->num_local_edges[0]; ++e) {
    std::set<int> p_l_colors;
    std::set<int> s_l_colors;
    for (int i=gs->eve_offsets[0][e]; i<gs->eve_offsets[0][e+1]; ++i) {
      p_l_colors.insert(gp->eve_lists[0][i]);
      s_l_colors.insert(gs->eve_lists[0][i]);
    }
    if (p_l_colors != s_l_colors)
      ++conflicts;
  }
  printf("Lists conflicts: %i\n", conflicts);

  destroyGraph(g_parallel);
  destroyGraph(g_serial);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}
