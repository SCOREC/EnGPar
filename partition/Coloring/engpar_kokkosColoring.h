#ifndef __ENGPAR_KOKKOS_COLORING__
#define __ENGPAR_KOKKOS_COLORING__

#include <agi.h>
#include <engpar_support.h>

namespace engpar {
  class ColoringInput;
  /** \brief return a host array with the coloring
   *  \remark this supports procedures that will use the
   *  coloring on the host
   */
  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t** colors);
}

#ifdef KOKKOS_ENABLED

#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosGraph_graph_color.hpp>
#include <KokkosKernels_Handle.hpp>

namespace engpar {
  /** \brief return a kokkos device view with the coloring
   *  \remark this supports procedures that will use the
   *  coloring on the device
   */
  LIDs EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t& numColors);
  agi::lid_t EnGPar_KokkosColoring(ColoringInput* in, agi::lid_t** colors);
}

#endif //KOKKOS_ENABLED
#endif // __ENGPAR_KOKKOS_COLORING__
