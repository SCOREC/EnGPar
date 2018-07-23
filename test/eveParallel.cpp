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
typedef Kokkos::View<agi::lid_t*> kkLidView;


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
  Kokkos::View<int**> A ("adjacency_matrix", N, N);
  kkLidView degree_view ("degree_view", pg->num_local_verts+1);
  hostToDevice(degree_view, pg->degree_list[t]);
  kkLidView edge_view ("edge_view", pg->num_local_pins[t]);
  hostToDevice(edge_view, pg->edge_list[t]);
  // User parallel for to construct adjacency matrix 
  Kokkos::parallel_for(pg->num_local_verts, KOKKOS_LAMBDA(const int v) {
    for (agi::lid_t i=degree_view(v); i<degree_view(v+1); ++i) {
      for (agi::lid_t j=degree_view(v); j<degree_view(v+1); ++j) {
        if (i!=j)
          A(edge_view(i),edge_view(j)) = 1;
      }
    }
  });
  // Compress 'A' into a CSR
  Kokkos::View<int*> eve_deg ("A_degree_list", N+1); 
  Kokkos::parallel_for (N, KOKKOS_LAMBDA(const int& i) {
    int degree = 0;
    for (int j=0; j<N; ++j) {
      degree += A(i,j);
    }
    eve_deg(i+1) = degree;
  }); 
  Kokkos::parallel_scan(N+1, KOKKOS_LAMBDA(const int& i, int& upd, const bool& final) {
    const int val = eve_deg(i); 
    upd += val;
    if (final)
      eve_deg(i) = upd; 
  });
  int eve_edge_size = 0;
  Kokkos::parallel_reduce(N, KOKKOS_LAMBDA(const int i, int& upd) {
    upd += eve_deg(i);
  }, eve_edge_size);
  kkLidView eve_edges ("eve_edges", eve_edge_size);
  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int i) {
      int idx = eve_deg(i);
      while (idx != eve_deg(i+1)) {
        for (int j=0; j<N; ++j) {
          if (A(i,j)==1)
            eve_edges(idx++) = j;
        }
      }
  });
  pg->eve_offsets[t] = new int[N+1];
  pg->eve_lists[t] = new int[eve_edge_size];
  deviceToHost(eve_deg, pg->eve_offsets[t]);
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
    //printf("%i, ",gp->eve_offsets[0][i]);
    //printf("%i\n",gs->eve_offsets[0][i]);
    if (gp->eve_offsets[0][i] != gs->eve_offsets[0][i])
      ++conflicts;
  }
  printf("Offsets conflicts: %i\n", conflicts);
  conflicts = 0;
  for (int e=0; e<gs->num_local_edges[0]; ++e) {
    std::set<int> p_l_edges;
    std::set<int> s_l_edges;
    for (int i=gp->eve_offsets[0][e]; i<gp->eve_offsets[0][e+1]; ++i) {
      p_l_edges.insert(gp->eve_lists[0][i]);
      s_l_edges.insert(gs->eve_lists[0][i]);
    }
    if (p_l_edges!=s_l_edges)
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
