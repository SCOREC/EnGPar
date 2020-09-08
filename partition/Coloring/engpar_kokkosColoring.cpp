#include <PCU.h>
#include <ngraph.h>
#include "engpar_kokkosColoring.h"
#include "engpar_coloring_input.h"

#ifdef KOKKOS_ENABLED
namespace engpar {
  LIDs EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t& numColors) { 
    agi::PNgraph* pg = in->g->publicize();
    // Retrieve relavent information from graph  
    agi::lid_t* adj_offsets = nullptr;
    agi::lid_t* adj_lists = nullptr;
    agi::lid_t numEnts = 0;
    if (in->primaryType == VTX_TYPE) {
      // Vertex colorings
      double t0 = PCU_Time();
      in->g->create_vev_adjacency(in->edgeType);
      numEnts = in->g->numLocalVtxs();
      adj_offsets = pg->vev_offsets[in->edgeType];
      adj_lists = pg->vev_lists[in->edgeType];
    }
    else {
      assert(in->g->isHyper());
      // Edge coloring
      double t0 = PCU_Time();
      //in->g->create_eve_adjacency(in->edgeType);
      in->g->parallel_create_eve(in->edgeType, in->boundaryOnly);
      numEnts = pg->num_local_edges[in->edgeType];
      adj_offsets = pg->eve_offsets[in->edgeType];
      adj_lists = pg->eve_lists[in->edgeType];
    }
    LIDs colors_d("colors_device", numEnts);
    // Create views
    LIDs adj_offsets_view ("adj_offsets_view", numEnts+1);
    hostToDevice(adj_offsets_view, adj_offsets);
    LIDs adj_lists_view ("adj_lists_view", adj_offsets[numEnts]);
    hostToDevice(adj_lists_view, adj_lists);
    // Typedefs to simplify kokkos template calls 
    typedef KokkosSparse::CrsMatrix<agi::lid_t, agi::lid_t, exeSpace::device_type, void, int> crsMat_t;
    typedef crsMat_t::StaticCrsGraphType graph_t;
    typedef graph_t::entries_type::non_const_type color_view_t;
    typedef graph_t::row_map_type lno_view_t;
    typedef graph_t::entries_type lno_nnz_view_t;
    typedef graph_t::entries_type::non_const_type  color_view_t;
    typedef KokkosKernels::Experimental::KokkosKernelsHandle <agi::lid_t, agi::lid_t, agi::lid_t, 
            exeSpace::execution_space, exeSpace::memory_space, exeSpace::memory_space> KernelHandle; 
    // Create kernel handle which will call graph color
    KernelHandle* kh = new KernelHandle();
    kh->set_team_work_size(32);
    kh->set_dynamic_scheduling(true);
    kh->create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
    // Run kokkos Coloring and delete handle
    double t0 = PCU_Time();
    KokkosGraph::Experimental::graph_color<KernelHandle, LIDs, LIDs>
      (kh, numEnts, numEnts, adj_offsets_view, adj_lists_view);
    colors_d = kh->get_graph_coloring_handle()->get_vertex_colors();
    numColors = kh->get_graph_coloring_handle()->get_num_colors();
    kh->destroy_graph_coloring_handle();  
    return colors_d;
  }

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t** colors) { 
    agi::lid_t numColors;
    LIDs colors_d = EnGPar_KokkosColoring(in, numColors);
    // Move coloring into array on host 
    *colors = new agi::lid_t[colors_d.dimension_0()];
    deviceToHost(colors_d, *colors);
    return numColors;
  }
}

#else

  agi::lid_t EnGPar_KokkosColoring(ColoringInput*, agi::lid_t*) {
    throw std::runtime_error("KOKKOS not found\n");
    return 0;
  }
}
#endif
