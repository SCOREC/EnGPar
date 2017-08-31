/*
//@HEADER
// *****************************************************************************
//
//       Multistep: (Strongly) Connected Components Algorithms
//              Copyright (2016) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  George M. Slota (gmslota@sandia.gov)
//                      Siva Rajamanickam (srajama@sandia.gov)
//                      Kamesh Madduri (madduri@cse.psu.edu)
//
// *****************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include "kokkosBfs.h"

#define out_vertice(n, j) out_array[out_degree_list[n]+j]
#define out_degree(n) (out_degree_list[n+1] - out_degree_list[n])

#if GPU
  #define LOCAL_BUFFER_LENGTH 16
  #define WORK_CHUNK 2048
  #define WARP_SIZE 32
#else
  #define LOCAL_BUFFER_LENGTH 2048
  #define WORK_CHUNK 2048
#endif
#define QUEUE_MULTIPLIER 1.25

typedef Kokkos::TeamPolicy<> team_policy;

typedef Kokkos::View<bool*> bool_array;
typedef Kokkos::View<int*> int_array;
typedef Kokkos::View<unsigned*> unsigned_array;
typedef Kokkos::View<int> int_type;
typedef Kokkos::View<double> double_type;

struct fwbw_init {
  int_type root;
  int_type n;
  bool_array fw;

  fwbw_init(int_type n_in,
    bool_array fw_in,
    int_type root_in) 
  : n(n_in)                   // device - # of vertices
  , fw(fw_in)                 // device - list of visited vertices
  , root(root_in)             // device - starting vertex
  {    
    typename int_type::HostMirror n_host = Kokkos::create_mirror(n);
    Kokkos::deep_copy(n_host, n);

    int num_teams = (n_host() + WORK_CHUNK - 1 ) / WORK_CHUNK;
    int team_size = team_policy::team_size_recommended(*this); //best guess
    team_policy policy(num_teams, team_size);
    Kokkos::parallel_for(policy , *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_policy::member_type &dev ) const
  {
    int begin = dev.league_rank() * WORK_CHUNK + dev.team_rank();
    int end = begin + WORK_CHUNK;
    end = n() < end ? n() : end;
    int team_size = dev.team_size();

    for (int v = begin; v < end; v += team_size)
    {
      if (v == root())
      {
        fw[v] = true;
      }
      else
      {
        fw[v] = false;
      }
    }
  }
};

struct fwbw_baseline {
  int_array out_array;
  int_array out_degree_list;
  bool_array valid;

  int num_visited;
  int_type root;

  int_array queue;
  int_array queue_next;
  int_type queue_size;
  int_type next_size;

  bool_array visited;

  fwbw_baseline(
    int_array out_degree_list_in,
    int_array out_array_in,
    bool_array valid_in,
    int_type root_in,
    bool_array visited_in,
    int_array queue_in,
    int_array queue_next_in)
  : out_array(out_array_in)             // device - list of adjacent vertices for each vtx
  , out_degree_list(out_degree_list_in) // device - degree of each vtx
  , valid(valid_in)                     // device - ignore, set to true 
  , root(root_in)                       // device - starting vtx
  , visited(visited_in)                 // device - list of visited vertices
  , queue(queue_in)                     // device - list of vertices to visit
  , queue_next(queue_next_in)           // device - list of vertices to visit next
  , queue_size("queue size")            // device - size of queue
  , next_size("next size")              // device - size of queue_next
  {
    typename int_type::HostMirror host_size = create_mirror(queue_size);
    typename int_type::HostMirror host_next = create_mirror(next_size);
    typename int_type::HostMirror host_root = create_mirror(root);
    deep_copy(host_root, root);
    host_size() = 1;
    num_visited = 1;
    deep_copy(queue_size, host_size);

    // initialize the queue with the root node
    typename int_array::HostMirror host_queue = create_mirror_view(queue);
    host_queue(0) = host_root();
    deep_copy(queue, host_queue);

    int team_size = team_policy::team_size_recommended(*this);
    int count = 0;
    while (host_size() > 0)
    {
      int num_teams = (host_size() + WORK_CHUNK - 1 ) / WORK_CHUNK;
      printf("iter %d num_visited %d host_size %d num_teams %d\n",
          count, num_visited, host_size(), num_teams);
      team_policy policy(num_teams, team_size);
      Kokkos::parallel_for(policy , *this);

      deep_copy(host_next, next_size);
      num_visited += host_next();
      host_size() = host_next();
      host_next() = 0;
      deep_copy(next_size, host_next);
      deep_copy(queue_size, host_size);

      int_array temp = queue;
      queue = queue_next;
      queue_next = temp;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_policy::member_type &dev ) const
  {
    int local_buffer[ LOCAL_BUFFER_LENGTH ];
    int local_count = 0;

    int begin = dev.league_rank() * WORK_CHUNK + dev.team_rank();
    int end = begin + WORK_CHUNK;
    int team_size = dev.team_size();

    for (int i = begin; i < end; i += team_size)
    {
      int out_degree = 0;
      int vert = -1;

      if (i < queue_size())
      {
        vert = queue[i];
        out_degree = out_degree(vert);
      }

      for (int j = 0; j < out_degree; ++j)
      {
        int out = out_vertice(vert, j);
        if (!visited[out] && valid[out])
        {
          visited[out] = true;
          local_buffer[local_count++] = out;
        }

        if (local_count == LOCAL_BUFFER_LENGTH)
        {
          const int thread_start = Kokkos::atomic_fetch_add(
            &next_size(), local_count);
          for (int l = 0; l < local_count; ++l)
            queue_next[thread_start+l] = local_buffer[l];
          local_count = 0;
        }
      }
    }

    int team_offset = dev.team_scan(local_count, &next_size());
    for (int l = 0; l < local_count; ++l)
      queue_next[team_offset+l] = local_buffer[l];
  }

};
