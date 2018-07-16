#include <iostream>
#include <PCU.h>
#include <ngraph.h>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosGraph_graph_color.hpp>
#include <KokkosKernels_Handle.hpp>
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
      // Edge coloring
      double t0 = PCU_Time();
      in->g->create_eve_adjacency(in->edgeType);
      printf ("eve partition time: %f\n", PCU_Time()-t0); 
      numEnts = in->g->numLocalEdges(in->edgeType);
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
    kh->set_team_work_size(16);
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

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t* colors) {
    throw std::runtime_error("KOKKOS not found\n");
    return 0;
  }

#endif
}
