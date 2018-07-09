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

  /** \brief helper function to transfer a device view to a host view
   */
  void deviceToHost(kkLidView d, kkLidView h) {
    Kokkos::deep_copy(h,d);
  }
  /** \brief helper function to transfer a device view to a host array
   */
  void deviceToHost(kkLidView d, agi::lid_t* h) {
    kkLidView::HostMirror hv = Kokkos::create_mirror_view(d);
    Kokkos::deep_copy(hv,d);
    for(size_t i=0; i<hv.size(); i++)
      h[i] = hv(i);
  }

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t** colors) {
    // FIXME - get graph entity count of in->primaryType
    agi::lid_t numEnts = in->g->numLocalVtxs();
    kkLidView colors_d("colors_device", numEnts);
    // FIXME - make this work with different edge types
    agi::PNgraph* pg = in->g->publicize();
    // Retrieve relavent information from graph
    const agi::lid_t numverts = pg->num_local_verts;
    const agi::lid_t numedges = pg->num_local_edges[0];
    const agi::lid_t* degree_list = pg->degree_list[0];
    const agi::lid_t* edge_list = pg->edge_list[0];
    // Create views
    kkLidView degree_view ("degree_view", numverts+1);
    kkLidView::HostMirror host_degree_view ("host_degree_view", numverts+1);
    kkLidView edge_view ("edge_view", numedges);
    kkLidView::HostMirror host_edge_view ("host_edge_view",numedges);
    for (int i=0; i<numverts+1; ++i) {
      host_degree_view(i) = degree_list[i];
    }
    for (int i=0; i<numedges; ++i) {
      host_edge_view(i) = edge_list[i];
    }
    Kokkos::deep_copy(degree_view, host_degree_view);
    Kokkos::deep_copy(edge_view, host_edge_view); 
    // Typedefs to simplify kokkos template calls
    typedef Kokkos::DefaultExecutionSpace exe_space;
    typedef KokkosSparse::CrsMatrix<agi::lid_t, agi::lid_t, exe_space::device_type, void, int> crsMat_t;
    typedef crsMat_t::StaticCrsGraphType graph_t;
    typedef graph_t::entries_type::non_const_type color_view_t;
    typedef graph_t::row_map_type lno_view_t;
    typedef graph_t::entries_type lno_nnz_view_t;
    typedef graph_t::entries_type::non_const_type  color_view_t;
    typedef KokkosKernels::Experimental::KokkosKernelsHandle <agi::lid_t, agi::lid_t, agi::lid_t, 
            exe_space, exe_space::memory_space, exe_space::memory_space> KernelHandle; 

    KernelHandle* kh = new KernelHandle();
    kh->set_team_work_size(16);
    kh->set_dynamic_scheduling(true);
    kh->create_graph_coloring_handle(KokkosGraph::COLORING_DEFAULT);
    // Run kokkos Coloring
    KokkosGraph::Experimental::graph_color<KernelHandle, kkLidView, kkLidView>
      (kh, numverts, numverts, degree_view, edge_view);
    colors_d = kh->get_graph_coloring_handle()->get_vertex_colors();
    kh->destroy_graph_coloring_handle(); 
    // Move coloring into array on host 
    *colors = new agi::lid_t[numEnts];
    deviceToHost(colors_d, *colors);
    return 0;
  }

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, kkLidView colors) {
    return 0;
  }

#else

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t* colors) {
    throw std::runtime_error("KOKKOS not found\n");
    return 0;
  }

#endif
}
