#ifndef NGRAPH_H__
#define NGRAPH_H__
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <map>
#include <cassert>
#include "agi.h"
#include "EdgeIterator.h"
#include "GraphIterator.h"

/** \file ngraph.h
    \brief The N-Graph interface */
namespace agi {

class GraphVertex;
class GraphEdge;
class VertexIterator;
class PinIterator;
/** \class Ngraph
    \brief An abstract graph used to represent the data passed into EnGPar
 
    Extending this class allows the user to fill the data structures with the user's data. Both regular graph and hypergraph designs are supported.
*/
class Ngraph {
  
public:
  Ngraph();
  virtual ~Ngraph();

  //Global Part Information
  /** \brief Returns the number of vertices in the graph across all processes */
  gid_t numGlobalVtxs() const {return num_global_verts;}
  /** \brief Returns the number of edges in the graph across all processors
   * \param t the edge type index */
  gid_t numGlobalEdges(etype t=0) const {return num_global_edges[t];}
  /** \brief Returns the number of pins in the graph across all processors [HG only]
   * \param t the edge type index */
  gid_t numGlobalPins(etype t=0) const {return num_global_pins[t];}
  
  //Local Part Information
  /** \brief Returns the number of vertices on this process */
  lid_t numLocalVtxs() const {return num_local_verts;}
  /** \brief Returns the number of ghost vertices on this process */
  lid_t numGhostVtxs() const {return num_ghost_verts;}
  /** \brief Returns the number of vertices+ghost vertices on this process */
  lid_t numTotalVtxs() const {return num_local_verts+num_ghost_verts;}
  /** \brief Returns the number of edge types in the graph */
  int numEdgeTypes() const {return num_types;}
  /** \brief Returns the number of edges on this process 
   * \param t the edge type index */
  lid_t numLocalEdges(etype t=0) const {return num_local_edges[t];}
  /** \brief Returns the number of pins on this process [HG only]
   * \param t the edge type index */
  lid_t numLocalPins(etype t=0) const {return num_local_pins[t];}
  /** \brief Returns whether the graph is in hypergraph mode or not*/
  bool isHyper() const {return isHyperGraph;}
  
  //Vertex Operations
  /** \brief Returns the weight of a vertex 
   * \param vtx the graph vertex
   * \return the weight
   */
  const wgt_t& weight(GraphVertex* vtx) const;
  /** \brief Returns the coordinates of a vertex
   * \param vtx the graph vertex
   * \return the coordinate
   */
  const coord_t& coord(GraphVertex* vtx) const;
  /** \brief Returns the owner of a vertex
   * \param vtx the graph vertex
   * \return the part id of the owner
   */
  part_t owner(GraphVertex* vtx) const;
  // \cond HACK
  lid_t localID(GraphVertex*) const;
  gid_t globalID(GraphVertex*) const;
  GraphVertex* find(GraphVertex* vtx) const;
  // \endcond
  
  //Edge Operations
    /** \brief Returns the weight of a edge 
   * \param edge the graph edge
   * \return the weight
   */
  double weight(GraphEdge* edge) const;
  // \cond HACK
  lid_t u(lid_t,etype t =0) const;
  // \endcond
  /** \brief Returns the source of an edge [Not HG]
   * \param edge the graph edge
   * \return the source vertex
   */
  GraphVertex* u(GraphEdge* edge) const;
  /** \brief Returns the destination of an edge [Not HG]
   * \param edge the graph edge
   * \return the destination vertex
   */
  GraphVertex* v(GraphEdge* edge) const;
  
  //Adjacency Operations
  /** \brief Retuns the degree of a vertex
   * \param vtx the graph vertex
   * \param t the edge type
   * \return the degree of the vertex
   */
  lid_t degree(GraphVertex* vtx,etype t=0) const;
  /** \brief Creates an iterator over the edges of the given vertex
   * \param vtx the graph vertex
   * \param t the edge type
   * \return an iterator that can loop over each edge of the vertex
   */
  EdgeIterator* edges(GraphVertex* vtx,etype t=0) const;
  /** \brief Create an iterator over neighboring vertices given the vertex
   *
   * \param vtx the graph vertex
   * \param t the edge type
   * \return an iterator that can loop over each neighboring vertex
   */
  GraphIterator* adjacent(GraphVertex* vtx, etype t=0) const;

  /** \brief Returns the degree of an edge [HG only]
   * \param edge the graph edge
   * \return the number of pins to vertices from this edge
   */
  lid_t degree(GraphEdge* edge) const;
  /** \brief Creates an iterator over the pins of the given hyperedge.
   * \param edge the graph hyperedge
   * \return an iterator that can loop over each vertex connected to this hyperedge
   */
  PinIterator* pins(GraphEdge* edge) const;
  
  //Iterator Traversal
  /** \brief Creates an iterator over all vertices
   * \return an iterator over all vertices in the graph
   */
  VertexIterator* begin() const;
  // \cond HACK
  GraphVertex* findGID(gid_t gid) const;
  // \endcond
  /** \brief Creates an iterator over all edges of a given type
   * \param t the edge type
   * \return an iterator over all edges of the given type
   */
  EdgeIterator* begin(etype t) const;
  /** \brief Iterates the vertex iterator
   * \param vitr a vertex iterator
   * \return the current graph vertex
   */
  GraphVertex* iterate(VertexIterator*& vitr) const;
  /** \brief Iterates the edge iterator
   * \param eitr the edge iterator
   * \return the current graph edge
   */
  GraphEdge* iterate(EdgeIterator*& eitr) const;
  /** \brief Iterates the pin iterator
   * \param pitr the pin iterator
   * \return the graph vertex at the end of the pin
   */
  GraphVertex* iterate(PinIterator*& pitr) const;
  /** \brief Iterates the graph iterator
   * \param gitr the graph iterator
   * \return the next adjacent graph vertex
   */
  GraphVertex* iterate(GraphIterator*& gitr) const;
  /** \brief Returns the current edge of the graph iterator
   * \param gitr the graph iterator
   * \return the edge
   */
  GraphEdge* edge(GraphIterator* gitr) const;
  //Destroys iterator
  /** \brief Cleans up the memory of an edge iterator
   * \param eitr the edge iterator
   */
  void destroy(EdgeIterator* eitr) const;
  /** \brief Cleans up the memory of a graph iterator
   * \param gitr the graph iterator
   */
  void destroy(GraphIterator* gitr) const;

  //Utility
  /** \brief Tests if two vertices are equal
     \param vtx1 the first vertex
     \param vtx2 the second vertex
  */
  bool isEqual(GraphVertex* vtx1,GraphVertex* vtx2) const;
  /** \brief Tests if two edges are equal
     \param edge1 the first edge
     \param edge2 the second edge
  */
  bool isEqual(GraphEdge* edge1,GraphEdge* edge2) const;
  /** \brief A method to provide a migration plan of the vertices
   * \param plan a map from graph vertex to part id
   */
  virtual void migrate(std::map<GraphVertex*,int>& plan) = 0;

 protected:
  /** \brief Adds an edge type to the graph
   * \return the new edge type id
   */
  etype addEdgeType() {assert(num_types!=MAX_TYPES);return num_types++;}

  // \cond
  void makeEdgeArray(etype t,int count);
  void setEdge(lid_t,gid_t,wgt_t,etype);
  // \endcond
  
  /*
  void create_csr(int num_verts, int num_edges, int* srcs,
                  int* dsts, int* wgts);
  */
  // \cond members
  //TODO: Try to compress global id

  bool isHyperGraph;
  //Global number of vertices and edges
  gid_t num_global_verts;
  gid_t num_global_edges[MAX_TYPES];
  gid_t num_global_pins[MAX_TYPES];
  
  // number of vertices and edges
  lid_t num_local_verts;
  lid_t num_ghost_verts;
  int num_types;
  lid_t num_local_edges[MAX_TYPES];
  lid_t num_local_pins[MAX_TYPES];
  
  // vertex weights
  // size = num_local_verts
  wgt_t* local_weights;
  // vertex coordinates
  // size = num_local_verts
  coord_t* local_coords;
  
  //edge weights
  // size=num_edges
  gid_t* edge_ids[MAX_TYPES];
  wgt_t* edge_weights[MAX_TYPES];
  
  
  lid_t* degree_list[MAX_TYPES];
  lid_t* edge_list[MAX_TYPES];
  lid_t* pin_degree_list[MAX_TYPES];
  lid_t* pin_list[MAX_TYPES];

  typedef std::unordered_map<gid_t,lid_t> map_t;
  map_t vtx_mapping;
  map_t edge_mapping[MAX_TYPES];
  //TODO: Tack ghost unmap on top of local_unmap
  gid_t* local_unmap;
  gid_t* ghost_unmap;
  part_t* owners;
  // \endcond
};
/** \brief Cleans up the memory of the graph
 * \param g the graph
 */
void destroyGraph(Ngraph* g);
} //namespace

#endif
