#include "engpar_reduce.h"
namespace engpar {
  typedef ENGPAR_LID_T lid_t;
  lid_t getMax(LIDs a) {
    lid_t max;
    const int n = a.extent(0);
    Kokkos::parallel_reduce(n, KOKKOS_LAMBDA(const int i, lid_t& lmax) {
      const lid_t val = a(i);
      if(val>lmax) {
        lmax = val;
      }
    },Kokkos::Max<lid_t>(max));
    return max;
  }
}
