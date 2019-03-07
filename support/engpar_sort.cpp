/** the following heavily borrows from 
 * Dan Ibanez's Omega_h sorting functions in
 * github.com/SNLComputation/omega_h.git
 * Omega_h_sort.[hpp|cpp]
 */

#include <Kokkos_Core.hpp>
#include <engpar_sort.h>
#include <algorithm>
#include <vector>

#if defined(ENGPAR_USE_CUDA)
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#elif defined(ENGPAR_USE_OPENMP)

#include <omp.h>
#include <pss/parallel_stable_sort.hpp>
#include <pss/pss_common.hpp>

#endif

namespace engpar {

template <typename T, typename Comp>
static void parallel_sort(T* b, T* e, Comp c) {
#if defined(ENGPAR_USE_CUDA)
  auto bptr = thrust::device_ptr<T>(b);
  auto eptr = thrust::device_ptr<T>(e);
  thrust::stable_sort(bptr, eptr, c);
#elif defined(ENGPAR_USE_OPENMP)
  pss::parallel_stable_sort(b, e, c);
#else
  std::stable_sort(b, e, c);
#endif
}

template <typename T>
struct CompareKeySets {
  T const* keys_;
  CompareKeySets(T const* keys) : keys_(keys) {}
  KOKKOS_INLINE_FUNCTION
  bool operator()(const int& a, const int& b) const {
    T x = keys_[a];
    T y = keys_[b];
    if (x != y) return x > y;
    return false;
  }
};

LIDs sort_by_keys(LIDs keys) {
  size_t n = keys.size();
  LIDs perm("permutation", n);
  Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int i) {
      perm(i) = i;
  });
  ENGPAR_LID_T* begin = perm.data();
  ENGPAR_LID_T* end = perm.data() + n;
  ENGPAR_LID_T const* keyptr = keys.data();
  typedef struct CompareKeySets<ENGPAR_LID_T> CompLid;
  CompLid c(keyptr);
  Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int i) {
      bool res = c(i,i+1);
      printf("%d %d %d\n", keyptr[i], keyptr[i+1], res);
  });
  parallel_sort<ENGPAR_LID_T, CompLid >(begin, end, c);
  return perm;
}

}
