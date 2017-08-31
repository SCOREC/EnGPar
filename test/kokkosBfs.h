#ifndef KOKKOSBFS_H
#define KOKKOSBFS_H

typedef Kokkos::TeamPolicy<> team_policy;

typedef Kokkos::View<bool*> bool_array;
typedef Kokkos::View<int*> int_array;
typedef Kokkos::View<unsigned*> unsigned_array;
typedef Kokkos::View<int> int_type;
typedef Kokkos::View<double> double_type;

#define QUEUE_MULTIPLIER 1.25
#define WORK_CHUNK 2048
#if GPU
  #define LOCAL_BUFFER_LENGTH 16
#else
  #define LOCAL_BUFFER_LENGTH 2048
#endif

#endif
