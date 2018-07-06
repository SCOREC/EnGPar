#include <iostream>
#include <PCU.h>
#include <ngraph.h>
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
    for(size_t i=0; i<hv.extent(0); i++)
      h[i] = hv(i);
  }

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t** colors) {
    agi::PNgraph* pg = in->g->publicize();
    agi::lid_t numEnts = 1; //FIXME - get graph entity count of in->primaryTyp
    kkLidView colors_d("colors_device", numEnts);
    agi::lid_t ret = EnGPar_KokkosColoring(in, colors_d);
    assert(ret);
    *colors = new agi::lid_t[numEnts];
    deviceToHost(colors_d, *colors);
    return 0;
  }

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, kkLidView colors) {
    // call kokkos coloring here and return the device view via the 'colors' arg
    return 0;
  }

#else

  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t* colors) {
    throw std::runtime_error("KOKKOS not found\n");
    return 0;
  }

#endif
}
