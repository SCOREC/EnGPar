#include "engpar_kokkosColoring.h"
#include <iostream>
#include <PCU.h>
#include <ngraph.h>
namespace engpar {
#ifdef KOKKOS_ENABLED

  /** \brief helper function to transfer a device view to a host array
   */
  kkLidView deviceToHost(kkLidView d) {
    return kkLidView(1);
  }

  agi::lid_t* EnGPar_KokkosColoring(ColoringInput* in) {
    kkLidView colors_d = EnGPar_KokkosColoring(in);
    kkLidView colors_h = deviceToHost(colors_d);
    // insert code here to copy a host view to an array
    return NULL;
  }

  kkLidView EnGPar_KokkosColoring(ColoringInput* in) {
    return kkLidView(1);
  }

#else

  agi::lid_t* EnGPar_KokkosColoring(ColoringInput* in) {
    throw std::runtime_error("KOKKOS not found\n");
    return NULL;
  }

#endif
}
