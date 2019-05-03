#include <PCU.h>
#include <ngraph.h>
#include "engpar_kokkosColoring.h"
#include "engpar_coloring_input.h"
namespace engpar {
#ifdef KOKKOS_ENABLED

  /** \brief helper function to transfer a host array to a device view
   */
  void hostToDevice(kkLidView d, agi::lid_t* h) {
    kkLidView::HostMirror hv = Kokkos::create_mirror_view(d);
    for (size_t i=0; i<hv.size(); ++i)
      hv(i) = h[i];
    Kokkos::deep_copy(d,hv);
  }
  /** \brief helper function to transfer a device view to a host array
   */
  void deviceToHost(kkLidView d, agi::lid_t* h) {
    kkLidView::HostMirror hv = Kokkos::create_mirror_view(d);
    Kokkos::deep_copy(hv,d); 
    for(size_t i=0; i<hv.size(); ++i)
      h[i] = hv(i);
  }


  void parallel_create_eve(agi::Ngraph* g, agi::etype t=0) {
    assert(g->isHyper());
    agi::PNgraph* pg = g->publicize();
    if (pg->eve_offsets[t]) {
      delete [] pg->eve_offsets[t];
      delete [] pg->eve_lists[t];
    }
    // Load graph info to device
    const int N = pg->num_local_edges[t];
    const int M = pg->num_local_verts;
    kkLidView degree_view ("degree_view", M+1);
    kkLidView edge_view ("edge_view", pg->num_local_pins[t]);
    hostToDevice(degree_view, pg->degree_list[t]); 
    hostToDevice(edge_view, pg->edge_list[t]); 
    // make hint for map size to avoid resize
    int numAdj = 0;
    Kokkos::parallel_reduce (M, KOKKOS_LAMBDA(const int v, int& upd) {
      for (int i=degree_view(v); i<degree_view(v+1); ++i) {
        for (int j=degree_view(v); j<degree_view(v+1); ++j) {
          if (edge_view(i)!=edge_view(j))
            upd++;
        }
      } 
    }, numAdj);
    // fill map
    Kokkos::UnorderedMap<Kokkos::pair<int,int>,void> m (numAdj);
    Kokkos::parallel_for (M, KOKKOS_LAMBDA(const int v) {
      for (int i=degree_view(v); i<degree_view(v+1); ++i) {
        for (int j=degree_view(v); j<degree_view(v+1); ++j) {
          if (edge_view(i)!=edge_view(j)) {
            m.insert( Kokkos::pair<int,int>(edge_view(i),edge_view(j)) );
          }
        }
      }
    });
    // create offset array
    kkLidView deg ("deg", N+1);
    Kokkos::parallel_for (m.hash_capacity(), KOKKOS_LAMBDA(const int i) {
      if ( m.valid_at(i) ) {
        Kokkos::pair<int,int> p = m.key_at(i);
        Kokkos::atomic_fetch_add( &deg(p.first+1), 1);
      }
    });
    Kokkos::parallel_scan (N, KOKKOS_LAMBDA (const int i, int& upd, const bool& final) {
      const int val = deg(i+1);
      upd += val;
      if (final)
        deg(i+1) = upd;
    });
    pg->eve_offsets[t] = new int[N+1];
    deviceToHost(deg, pg->eve_offsets[t]);
    // make edge list array
    kkLidView edgeList ("edgeList", pg->eve_offsets[t][N]);
    kkLidView adjCount ("adjCount", N);
    Kokkos::parallel_for (m.hash_capacity(), KOKKOS_LAMBDA(const int i) {
      if ( m.valid_at(i) ) {
        Kokkos::pair<int,int> p = m.key_at(i);
        int e = deg(p.first);
        int idx = Kokkos::atomic_fetch_add( &adjCount(p.first), 1);
        edgeList(e+idx) = p.second;
      }
    });
    pg->eve_lists[t] = new int[pg->eve_offsets[t][N]];
    deviceToHost(edgeList, pg->eve_lists[t]);
  }


  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t** colors) { 
    agi::PNgraph* pg = in->g->publicize();
    // Retrieve relavent information from graph  
    agi::lid_t* adj_offsets = nullptr;
    agi::lid_t* adj_lists = nullptr;
    agi::lid_t numEnts = 0;
    if (in->primaryType == VTX_TYPE) {
      // Vertex colorings
      double t0 = PCU_Time();
      in->g->create_vev_adjacency(in->edgeType);
      printf ("eve partition time: %f\n", PCU_Time()-t0); 
      numEnts = in->g->numLocalVtxs();
      adj_offsets = pg->vev_offsets[in->edgeType];
      adj_lists = pg->vev_lists[in->edgeType];
    }
    else {
      assert(in->g->isHyper());
      // Edge coloring
      double t0 = PCU_Time();
      //in->g->create_eve_adjacency(in->edgeType);
      parallel_create_eve(in->g, in->edgeType);
      printf ("eve partition time: %f\n", PCU_Time()-t0); 
      numEnts = pg->num_local_edges[in->edgeType];
      adj_offsets = pg->eve_offsets[in->edgeType];
      adj_lists = pg->eve_lists[in->edgeType];
    }
    kkLidView colors_d("colors_device", numEnts);
    // Create views
    kkLidView adj_offsets_view ("adj_offsets_view", numEnts+1);
    hostToDevice(adj_offsets_view, adj_offsets);
    kkLidView adj_lists_view ("adj_lists_view", adj_offsets[numEnts]);
    hostToDevice(adj_lists_view, adj_lists);
    // Typedefs to simplify kokkos template calls 
    typedef KokkosSparse::CrsMatrix<agi::lid_t, agi::lid_t, exe_space::device_type, void, int> crsMat_t;
    typedef crsMat_t::StaticCrsGraphType graph_t;
    typedef graph_t::entries_type::non_const_type color_view_t;
    typedef graph_t::row_map_type lno_view_t;
    typedef graph_t::entries_type lno_nnz_view_t;
    typedef graph_t::entries_type::non_const_type  color_view_t;
    typedef KokkosKernels::Experimental::KokkosKernelsHandle <agi::lid_t, agi::lid_t, agi::lid_t, 
            exe_space::execution_space, exe_space::memory_space, exe_space::memory_space> KernelHandle; 
    // Create kernel handle which will call graph color
    KernelHandle* kh = new KernelHandle();
    kh->set_team_work_size(32);
    kh->set_dynamic_scheduling(true);
    kh->create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
    // Run kokkos Coloring and delete handle
    double t0 = PCU_Time();
    KokkosGraph::Experimental::graph_color<KernelHandle, kkLidView, kkLidView>
      (kh, numEnts, numEnts, adj_offsets_view, adj_lists_view);
    printf ("Coloring time: %f\n", PCU_Time()-t0);
    colors_d = kh->get_graph_coloring_handle()->get_vertex_colors();
    kh->destroy_graph_coloring_handle();  
    // Move coloring into array on host 
    *colors = new agi::lid_t[numEnts];
    deviceToHost(colors_d, *colors);
    return 0;
  }
 
#else

  agi::lid_t EnGPar_KokkosColoring(ColoringInput*, agi::lid_t*) {
    throw std::runtime_error("KOKKOS not found\n");
    return 0;
  }

#endif
}
