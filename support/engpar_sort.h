/** the following heavily borrows from 
 * Dan Ibanez's Omega_h sorting functions in
 * github.com/SNLComputation/omega_h.git
 * Omega_h_sort.[hpp|cpp]
 */

#ifndef ENGPAR_SORT_H
#define ENGPAR_SORT_H
#include <engpar_support.h>

/* Compute the permutation which sorts the given keys.
   Each key is a tuple of N integers of type T.
   T may be 32 or 64 bits, and N may be 1, 2, or 3.
   Let perm = sort_by_keys(keys, width);
   Then the key at perm[i] is <= the key at perm[i + 1].
   In other words, old_key_index = perm[new_key_index]
   The sorting algorithm used is stable, so if two keys
   are equal then their relative order after permutation
   will be the same as before the permutation.
   Tuples (keys) are sorted into lexical order, so they
   will be sorted by the first integer first, second
   integer second, etc.
 */
namespace engpar {
LIDs sort_by_keys(LIDs keys);
}

#endif
