#ifndef NGRAPH_H__
#define NGRAPH_H__
#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <cassert>
#include "agi.h"


/** \file ngraph.h
    \brief The N-Graph interface */
namespace agi {

class GraphVertex;
class GraphEdge;
class VertexIterator;
class PinIterator;
class EdgeIterator;
class GraphIterator;

/** \class Ngraph
    \brief An abstract graph used to represent the data passed into EnGPar
 
    Extending this class allows the user to fill the data structures with the user's data. Both regular graph and hypergraph designs are supported.
*/
class Ngraph {
  
public:
  Ngraph();
  /** \brief Constructs the Ngraph given a set of information
   * \param isHG true if the given construction is for a hypergraph
   * \param verts list of global ids of vertices that this part owns
   * \param weights list of the weights of each vertex
   * \param edge_ids list of global ids of edges that this part has
   * \param degs list of degrees of each edge (always 2 if 
   *        constructing a traiditional graph)
   * \param pins_to_verts list of the vertices the edges are connected to
   * \param owns mapping from global_id to owner for each ghosted vertex
   */
  void constructGraph(bool isHG,
		      std::vector<gid_t>& verts,
		      std::vector<wgt_t>& weights,
		      std::vector<gid_t>& edge_ids,
		      std::vector<lid_t>& degs,
		      std::vector<gid_t>& pins_to_verts,
		      std::unordered_map<gid_t,part_t>& owns);
  /** \brief Constructs the vertices of the Ngraph
   * \param isHG true if the given construction is for a hypergraph
   * \param verts list of global ids of vertices that this part owns
   * \param weights list of the weights of each vertex
   *
   * Must be called before constructEdges and should only be called once
   */
  void constructVerts(bool isHG,
		      std::vector<gid_t>& verts,
		      std::vector<wgt_t>& weights);
  /** \brief Constructs an edge type and returns the id of the type
   * \param edge_ids list of global ids of edges that this part has
   * \param degs list of degrees of each edge (always 2 if 
   *        constructing a traiditional graph)
   * \param pins_to_verts list of the vertices the edges are connected to
   *
   * Must be called after constructVerts and should be called once per edge type
   */
  etype constructEdges(std::vector<gid_t>& edge_ids,
		      std::vector<lid_t>& degs,
		      std::vector<gid_t>& pins_to_verts);
  /** \brief Constructs the ghost information for all non local vertices connected by edges
   * \param owns mapping from global_id to owner for each ghosted vertex
   *
   * Must be called after all edge types have been constructed 
   */
  void constructGhosts(std::unordered_map<gid_t,part_t>& owns);

  virtual ~Ngraph();
  // \cond
  void destroyData();
  // \endcond

  
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
  // \cond
  part_t originalOwner(GraphVertex* vtx) const;
  void setOriginalOwners();
  void setOriginalOwners(std::vector<part_t>&);
  // \endcond
  void getResidence(GraphEdge* e,Peers& residence) const;
  

  lid_t localID(GraphVertex*) const;
  gid_t globalID(GraphVertex*) const;
  // \cond HACK
  GraphVertex* find(GraphVertex* vtx) const;
  // \endcond
  
  //Edge Operations
    /** \brief Returns the weight of a edge 
   * \param edge the graph edge
   * \return the weight
   */
  double weight(GraphEdge* edge) const;
  
  lid_t localID(GraphEdge*) const;
  gid_t globalID(GraphEdge*) const;
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
  //Destroys iterator
  /** \brief Cleans up the memory of an edge iterator
   * \param eitr the edge iterator
   */
  void destroy(PinIterator* eitr) const;
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
  virtual PartitionMap* getPartition();
  // \cond
  /** \brief Sets the weights of the edges of a specific type
   * \param wgts the edge weights
   * \param t the edge type of the weights
   */
  void setEdgeWeights(std::vector<wgt_t>& wgts, etype t);

  void migrate(Migration* plan);
  // \endcond


  void saveToFile(char* prefix);
  void loadFromFile(char* prefix);
  
 protected:
  /** \brief Adds an edge type to the graph
   * \return the new edge type id
   */
  etype addEdgeType() {assert(num_types!=MAX_TYPES);return num_types++;}

  /** \brief Creates the edge array for a certain type
   * \param t the edge type
   * \param count the number of edges of this type
   */
  void makeEdgeArray(etype t,int count);
  /** \brief sets the values of a specific edge
   * \param l the local id of the edge
   * \param g the global id of the edge
   * \param w the weight of the edge
   * \param t the type of the edge
   */
  void setEdge(lid_t l,gid_t g,wgt_t w,etype t);

  /*
  void create_csr(int num_verts, int num_edges, int* srcs,
                  int* dsts, int* wgts);
  */
  // \cond DEV
  //TODO: Try to compress global id
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
  /** \brief A mapping from local id to global id for edges of each type
   *
   * size = num_local_edges
   */
  gid_t* edge_unmap[MAX_TYPES];  

  // \endcond
  // \cond
 private:
  void updateGhostOwners(Migration*);
  void sendVertex(GraphVertex*, part_t);
  void recvVertex(std::vector<gid_t>&,std::vector<wgt_t>&,
                  std::vector<part_t>&);
  void sendEdges(Migration*,std::unordered_set<GraphEdge*>&);
  void recvEdges(std::unordered_set<gid_t>&,std::vector<gid_t>&,
                 std::vector<wgt_t>&,std::vector<lid_t>&,
                 std::vector<gid_t>&,std::unordered_map<gid_t,part_t>&,
                 etype);
  // \endcond
};
/** \brief Cleans up the memory of the graph
 * \param g the graph
 */
void destroyGraph(Ngraph* g);
} //namespace


#ifdef KOKKOS_ENABLED
#include <Kokkos_Core.hpp>
#define KOKKOS_FOR_VERTS(g,v)						\
  Kokkos::parallel_for(g->numLocalVtxs(),KOKKOS_LAMBDA(uintptr_t i) {	\
      agi::GraphVertex* v = reinterpret_cast<agi::GraphVertex*>((char*)(i+1));
      
#define KOKKOS_END_FOR() });
#endif


#endif
