#ifndef __ENGPAR_KOKKOS_COLORING__
#define __ENGPAR_KOKKOS_COLORING__

namespace engpar {
  class ColoringInput;
  /** \brief return a host array with the coloring
   *  \remark this supports procedures that will use the
   *  coloring on the host
   */
  agi::lid_t* EnGPar_KokkosColoring(ColoringInput*);
}

#ifdef KOKKOS_ENABLED

#include <Kokkos_Core.hpp>
namespace engpar {
  typedef Kokkos::View<agi::lid_t*> kkLidView;

  /** \brief return a kokkos device view with the coloring
   *  \remark this supports procedures that will use the
   *  coloring on the device
   */
  kkLidView EnGPar_KokkosColoring(ColoringInput*);
}

#endif
