#ifndef __ENGPAR_KOKKOS_COLORING__
#define __ENGPAR_KOKKOS_COLORING__

#include <agi.h>

namespace engpar {
  class ColoringInput;
  /** \brief return a host array with the coloring
   *  \remark this supports procedures that will use the
   *  coloring on the host
   */
  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t** colors);
}

#ifdef KOKKOS_ENABLED

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosGraph_graph_color.hpp>
#include <KokkosKernels_Handle.hpp>

namespace engpar {
  typedef Kokkos::DefaultExecutionSpace exe_space;
  typedef Kokkos::View<agi::lid_t*, exe_space::device_type> kkLidView;
  typedef Kokkos::TeamPolicy<>::member_type member_type;

  /** \brief return a kokkos device view with the coloring
   *  \remark this supports procedures that will use the
   *  coloring on the device
   */
  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, kkLidView colors);
}

#endif //KOKKOS_ENABLED
#endif // __ENGPAR_KOKKOS_COLORING__
