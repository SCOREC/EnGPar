#ifndef __PUBLIC_NGRAPH_H__
#define __PUBLIC_NGRAPH_H__

#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <cassert>
#include "agi.h"

namespace agi {

  class GraphVertex;
  class GraphEdge;
  
class PNgraph {
 public:
    /** \brief A flag for if the hypergraph functionality is turned on
   */
  bool isHyperGraph;

  /** \brief The number of edge types
   */
  int num_types;

  /** \brief The number of unique vertices across all parts
   */
  gid_t num_global_verts;
  /** \brief The number of unique edges of each type across all parts
   */
  gid_t num_global_edges[MAX_TYPES];
  /** \brief The number of unique pins of each type across all parts
   * 
   * Only needed when isHyperGraph=true.
   */
  gid_t num_global_pins[MAX_TYPES];
  
  /** \brief The number of vertices on this part
   */
  lid_t num_local_verts;
  /** \brief The number of ghosted vertices on the part
   */
  lid_t num_ghost_verts;
  /** \brief The number of edges on the part of each type
   *
   * This includes edges going from local vertices to ghosted vertices
   */
  lid_t num_local_edges[MAX_TYPES];
  /** \brief The number of pins on the part of each type
   *
   * Only needed when isHyperGraph=true.
   * This includes pins going from edges to ghost vertices.
   */
  lid_t num_local_pins[MAX_TYPES];

  /** \brief The original owners of each vertex
   *
   * size = num_local_verts
   */
  part_t* original_owners;
  /** \brief The weights of the vertices.
   *
   * size = num_local_verts
   */
  wgt_t* local_weights;
  /** \brief The coordinates of the vertices
   *
   * size = num_local_verts
   * Only needed if using a geometric partitioner
   */
  coord_t* local_coords;
  
  /** \brief The weights of the edges
   *
   * size = num_local_edges
   */
  wgt_t* edge_weights[MAX_TYPES];
  
  /** \brief A stricly increasing list of the offsets into edge_list for each vertex for each edge type
   *
   * size = num_local_verts+1
   * Each entry represents the offset into the beginning of the vertex's edges.
   * The degree is calculated by: degree_list[type][i+1] - degree_list[type][i]
   */
  lid_t* degree_list[MAX_TYPES];
  /** \brief The list of edges off of each vertex
   *
   * size = num_local_edges
   * The starting position for a given vertex is found at degree_list[type][vertex]
   * and goes up to degree_list[type][vertex+1].
   */
  lid_t* edge_list[MAX_TYPES];
  /** \brief A stricly increasing list of the offsets into pin_list for each edge for each edge type
   *
   * size = num_local_edges+1
   * Only needed when isHyperGraph=true.
   * Each entry represents the offset into the beginning of the edges's pins.
   * The degree is calculated by: pin_degree_list[type][i+1] - pin_degree_list[type][i]
   */
  lid_t* pin_degree_list[MAX_TYPES];
  /** \brief The list of pins off of each edge
   *
   * size = num_local_pins
   * Only needed when isHyperGraph=true.
   * The starting position for a given edge is found at pin_degree_list[type][edge]
   * and goes up to pin_degree_list[type][edge+1].
   */

  lid_t* pin_list[MAX_TYPES];

  /** \brief A typedef for the mapping from global id to local id
   */
  typedef std::unordered_map<gid_t,lid_t> map_t;
  /** \brief A mapping from global id to local id for vertices
   *
   * Includes ghost vertices mapping
   */
  map_t vtx_mapping;
  /** \brief A mapping from local id to global id for vertices
   *
   * size = num_local_verts
   */
  gid_t* local_unmap;
  /** \brief A mapping from local id to global id for ghost vertices
   *
   * size = num_ghost_verts
   */
  gid_t* ghost_unmap;
  /** \brief The part owners of each ghost vertex
   *
   * size = num_ghost_verts
   */
  part_t* owners;

  /** \brief A mapping from global id to local id for edges of each type
   */

  map_t edge_mapping[MAX_TYPES];

  // \endcond
  // \cond DEV
  /** \brief A mapping from local id to global id for edges of each type
   *
   * size = num_local_edges
   */
  gid_t* edge_unmap[MAX_TYPES];  


  /* Second Adjacencies */

  /** \brief A stricly increasing list of the offsets for vertex-edge-vertex adjacencies
   *
   * size = num_local_verts+1
   * Each entry represents the offset into the beginning of the vertex's second adjacencies.
   * The 'degree' is calculated by: vev_list[type][i+1] - vev_list[type][i]
   */
  lid_t* vev_offsets[MAX_TYPES];
  /** \brief The list of vertex-edge-vertex adjacencies for each vertex
   *
   * size = vev_offsets[type][num_local_verts]
   * The starting position for a given vertex is found at vev_offsets[type][vertex]
   * and goes up to vev_offsets[type][vertex+1].
   */
  lid_t* vev_lists[MAX_TYPES];

  /** \brief A stricly increasing list of the offsets for edge-vertex-edge adjacencies
   *
   * size = num_local_edges[type]+1
   * Each entry represents the offset into the beginning of the edges's second adjacencies.
   * The 'degree' is calculated by: eve_list[type][i+1] - eve_list[type][i]
   */
  lid_t* eve_offsets[MAX_TYPES];
  /** \brief The list of edge-vertex-edge adjacencies for each edge
   *
   * size = eve_offsets[type][num_local_edges[type]]
   * The starting position for a given edge is found at eve_offsets[type][edge]
   * and goes up to eve_offsets[type][edge+1].
   */
  lid_t* eve_lists[MAX_TYPES];
  
  GraphVertex* getVertex(lid_t lid);
  GraphEdge* getEdge(lid_t lid,etype t);
};

}
#endif
