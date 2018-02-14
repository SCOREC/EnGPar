#ifndef NGRAPH_H__
#define NGRAPH_H__

#include <stdexcept>
#include "pngraph.h"
#include "agiMigration.h"
/** \file ngraph.h
    \brief The N-Graph interface 
*/
namespace agi {

class GraphVertex;
class GraphEdge;
class VertexIterator;
class GhostIterator;
class PinIterator;
class EdgeIterator;
class GraphIterator;
class GraphTag;

class MigrationTimers;

/** \class Ngraph
 *  \brief An abstract graph used to represent the data passed into EnGPar
 *   \nosubgrouping
 *
 *  Extending this class allows the user to fill the data structures with the user's data. Both regular graph and hypergraph designs are supported.
 *  Alternatively one can run the construct functions to build a graph without the use of an interface class.
*/
 class Ngraph : protected PNgraph {
  
public:
  
  /** \name Construction */
  ///@{
  friend Ngraph* createEmptyGraph();
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
  void constructVerts(bool isHG, lid_t num_verts,
                      gid_t* verts, wgt_t* weights = NULL);
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
                       std::vector<gid_t>& pins_to_verts,
                       std::vector<wgt_t>& e_weights);
  etype constructEdges(gid_t num_edges, gid_t* edge_ids,
                       lid_t* degs, gid_t* pins_to_verts,
                        wgt_t* e_weights = NULL);
  /** \brief Constructs the ghost information for all non local vertices connected by edges
   * \param owns mapping from global_id to owner for each ghosted vertex
   *
   * Must be called after all edge types have been constructed 
   */
  void constructGhosts(std::unordered_map<gid_t,part_t>& owns);
  void constructGhosts(lid_t num_ghosts, gid_t* vert_ids, part_t* owns);
  /** \brief Constructs the Ngraph given a set of information
   * \param isHG true if the given construction is for a hypergraph
   * \param verts list of global ids of vertices that this part owns
   * \param weights list of the weights of each vertex
   * \param edge_ids list of global ids of edges that this part has
   * \param degs list of degrees of each edge (always 2 if 
   *        constructing a traiditional graph)
   * \param pins_to_verts list of the vertices the edges are connected to
   * \param owns mapping from global_id to owner for each ghosted vertex
   *
   * This function is equivalent to calling constructVerts,constructEdges, and constructGhosts.
   */
  void constructGraph(bool isHG,
                      std::vector<gid_t>& verts,
                      std::vector<wgt_t>& weights,
                      std::vector<gid_t>& edge_ids,
                      std::vector<lid_t>& degs,
                      std::vector<gid_t>& pins_to_verts,
                      std::unordered_map<gid_t,part_t>& owns);
  /** \brief Removes all edges of type t from the graph.
   * \param t the edge type to be removed
   *
   * The function removes the memory of an edge type. This is useful to reduce the memory needed to transfer during migrations if the edges are no longer needed.
   */
  void removeEdges(etype t);
  
  // \cond
  PNgraph* publicize() {return this;}
  virtual ~Ngraph();
  // \endcond

  /** \brief Saves the graph connectivity information to prefix_#.bgd
   * \param prefix The prefix for the name of the file(s)
   *
   * One file is made per process
   */
  void saveToFile(const char* prefix);
  //TODO: put checks for number of files = number of processes
  /** \brief Loads the graph connectivity information from prefix_#.bgd
   * \param prefix The prefix of the name of the file(s)
   *
   * There should be one file per process
   */
  void loadFromFile(const char* prefix);
  ///@}


  /** \name Global Information */
  ///@{
  /** \brief Returns the number of vertices in the graph across all processes */
  gid_t numGlobalVtxs() const {return num_global_verts;}
  /** \brief Returns the number of edges in the graph across all processors
   * \param t the edge type index */
  gid_t numGlobalEdges(etype t=0) const {return num_global_edges[t];}
  /** \brief Returns the number of pins in the graph across all processors [HG only]
   * \param t the edge type index */
  gid_t numGlobalPins(etype t=0) const {return num_global_pins[t];}
  /** \brief Returns whether the graph is in hypergraph mode or not*/
  bool isHyper() const {return isHyperGraph;}
  ///@}
  /** \name Local Information */
  ///@{
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
  ///@}
  
  /** \name Vertex Operations */
  ///@{
  /** \brief Returns the weight of a vertex 
   * \param vtx the graph vertex
   * \return the weight
   */
  const wgt_t& weight(GraphVertex* vtx) const;
  /** \brief Returns true if coordinates are attached to the vertices.
   * \return whether coordinates are provided or not.
   */
  bool hasCoords() const;
  /** \brief Returns the coordinates of a vertex
   * \param vtx the graph vertex
   * \return the coordinate
   */
  const coord_t& coord(GraphVertex* vtx) const;
  /** \brief Set the coordinates of the vertices
   * \param cs an array of 3-D coordinates of size numLocalVtxs()
   */
  void setCoords(coord_t* cs);
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
  /** \brief Determines the parts that the given edge resides on
   * \param e the edge 
   * \param residence The part ids that the edge resides on
   */
  void getResidence(GraphEdge* e,Peers& residence) const;
  /** \brief Returns true if the target part shares the given edge
   * \param e the edge
   * \param peer the target part
   * \return whether the target part shares the edge
   */
  bool isResidentOn(GraphEdge* e,part_t peer) const;

  /** \brief Returns the local id of the vertex
   * \param vtx The vertex
   * \return The local id of the vertex.
   */
  lid_t localID(GraphVertex* vtx) const;
  /** \brief Returns the global id of the vertex.
   * \param vtx The vertex
   * \return The global id of the vertex.
   */
  gid_t globalID(GraphVertex* vtx) const;
  // \cond HACK
  GraphVertex* find(GraphVertex* vtx) const;
  // \endcond
  ///@}
  /** \name Edge Operations */
  ///@{
  /** \brief Returns the weight of a edge.
   * \param edge The graph edge
   * \return The weight
   */
  double weight(GraphEdge* edge) const;
  /** \brief Returns the local id of the edge.
   * \param e The edge
   * \return The local id of the edge.
   */  
  lid_t localID(GraphEdge* e) const;
  /** \brief Returns the global id of the edge.
   * \param e The edge
   * \return The global id of the edge.
   */  
  gid_t globalID(GraphEdge* e) const;
  // \cond
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
  /** \brief Sets the weights of the edges of a specific type
   * \param wgts the edge weights
   * \param t the edge type of the weights
   */
  void setEdgeWeights(std::vector<wgt_t>& wgts, etype t);
  ///@}
  
  /** \name Adjacency Operations */
  ///@{
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

  /** \brief Returns the degree of an edge
   * \param edge the graph edge
   * \return the number of pins to vertices from this edge
   */
  lid_t degree(GraphEdge* edge) const;
  /** \brief Creates an iterator over the pins of the given hyperedge.
   * \param edge the graph hyperedge
   * \return an iterator that can loop over each vertex connected to this hyperedge
   */
  PinIterator* pins(GraphEdge* edge) const;
  ///@}

  /** \name Tag Data */
  ///@{
  /** \brief Creates and returns a tag of integers over a type of entity
   * \param t the entity type (nothing defaults to vertices)
   * \return The tag to the data
   */
  GraphTag* createIntTag(etype t=VTX_TYPE);
  /** \brief Creates and returns a tag of doubles over a type of entity
   * \param t the entity type (nothing defaults to vertices)
   * \return The tag to the data
   */
  GraphTag* createDoubleTag(etype t=VTX_TYPE);
  /** \brief Creates and returns a tag of longs over a type of entity
   * \param t the entity type (nothing defaults to vertices)
   * \return The tag to the data
   */
  GraphTag* createLongTag(etype t=VTX_TYPE);
  /** \brief Destroys a tag created earlier
   * /param t The tag to be deleted
   */
  void destroyTag(GraphTag* t);
  /** \brief Retrieves an integer value of a tag for a vertex.
   * \param tag The tag.
   * \param v The vertex.
   * \return The value of the vertex for this tag 
   */
  int getIntTag(GraphTag* tag,GraphVertex* v);
  /** \brief Retrieves an integer value of a tag for a edge.
   * \param tag The tag.
   * \param e The edge.
   * \return The value of the edge for this tag 
   */
  int getIntTag(GraphTag* tag,GraphEdge* e);
  /** \brief Retrieves a double value of a tag for a vertex.
   * \param tag The tag.
   * \param v The vertex.
   * \return The value of the vertex for this tag 
   */
  double getDoubleTag(GraphTag* tag,GraphVertex* v);
  /** \brief Retrieves a double value of a tag for a edge.
   * \param tag The tag.
   * \param e The edge.
   * \return The value of the edge for this tag 
   */
  double getDoubleTag(GraphTag* tag,GraphEdge* e);
  /** \brief Retrieves a long value of a tag for a vertex.
   * \param tag The tag.
   * \param v The vertex.
   * \return The value of the vertex for this tag 
   */
  long getLongTag(GraphTag* tag,GraphVertex* v);
  /** \brief Retrieves a long value of a tag for a edge.
   * \param tag The tag.
   * \param e The edge.
   * \return The value of the edge for this tag 
   */
  long getLongTag(GraphTag* tag,GraphEdge* e);
  /** \brief Sets an integer value of a tag for a vertex.
   * \param tag The tag.
   * \param v The vertex.
   * \param val The value to be assigned to the vertex for the given tag
   */
  void setIntTag(GraphTag* tag,GraphVertex* v,int val);
  /** \brief Sets an integer value of a tag for a edge.
   * \param tag The tag.
   * \param e The edge.
   * \param val The value to be assigned to the edge for the given tag
   */
  void setIntTag(GraphTag* tag,GraphEdge* e,int val);
  /** \brief Sets a double value of a tag for a vertex.
   * \param tag The tag.
   * \param v The vertex.
   * \param val The value to be assigned to the vertex for the given tag
   */
  void setDoubleTag(GraphTag* tag,GraphVertex* v,double val);
  /** \brief Sets a double value of a tag for a edge.
   * \param tag The tag.
   * \param e The edge.
   * \param val The value to be assigned to the edge for the given tag
   */
  void setDoubleTag(GraphTag* tag,GraphEdge* e,double val);
  /** \brief Sets a long value of a tag for a vertex.
   * \param tag The tag.
   * \param v The vertex.
   * \param val The value to be assigned to the vertex for the given tag
   */
  void setLongTag(GraphTag* tag,GraphVertex* v,long val);
  /** \brief Sets a long value of a tag for a edge.
   * \param tag The tag.
   * \param e The edge.
   * \param val The value to be assigned to the edge for the given tag
   */
  void setLongTag(GraphTag* tag,GraphEdge* e,long val);
  ///@}
  
  /** \name Iterator Traversal*/
  ///@{
  /** \brief Creates an iterator over local vertices
   * \return an iterator over all vertices in the graph
   */
  VertexIterator* begin() const;
  /** \brief Creates an iterator over ghost vertices
   * \return an iterator over ghost vertices in the graph
   */  
  GhostIterator* beginGhosts() const;
  // \cond 
  GraphVertex* findGID(gid_t gid) const;
  // \endcond
  /** \brief Creates an iterator over all edges of a given type
   * \param t the edge type
   * \return an iterator over all edges of the given type
   */
  EdgeIterator* begin(etype t) const;
  /** \brief Iterates the vertex iterator
   * \param vitr a vertex iterator
   * \return the current vertex
   */
  GraphVertex* iterate(VertexIterator*& vitr) const;
  /** \brief Iterates the ghost vertex iterator
   * \param vitr a ghost vertex iterator
   * \return the current ghost vertex
   */
  GraphVertex* iterate(GhostIterator*& vitr) const;

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
  ///@}
  /** \name Utility */
  ///@{
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
  ///@}

  /** \name Partition/Migration */
  ///@{
  /** \brief Retrieves a migration plan for use on the original data
   * \return a map from global id to part id
   */
  virtual PartitionMap* getPartition();

  // \cond
  void migrate(Migration* plan, MigrationTimers* mt = NULL);
  // \endcond
  void repartition(part_t* partition);
  ///@}
  void changeOwners(int* newRanks);
  

 protected:
  // \cond 
  Ngraph();
  Ngraph(const Ngraph&) {throw std::runtime_error("Copying Ngraph not supported");}
  Ngraph& operator=(const Ngraph&) {
    throw std::runtime_error("Assignment of Ngraph not supported");
  }
  void destroyData();
  // \endcond

  // \cond  
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

  // \endcond
  
  // \cond PRIVATE
 private:
  void updateGhostOwners(Migration*);
  void sendVertex(GraphVertex*, part_t);
  void recvVertex(std::vector<gid_t>&,std::vector<wgt_t>&,
                  std::vector<part_t>&);
  void sendCoord(GraphVertex*, part_t);
  void recvCoord(coord_t*,lid_t& size);
  void sendEdges(Migration*,std::unordered_set<GraphEdge*>&);
  void recvEdges(std::unordered_set<gid_t>&,std::vector<gid_t>&,
                 std::vector<wgt_t>&,std::vector<lid_t>&,
                 std::vector<gid_t>&,std::unordered_map<gid_t,part_t>&,
                 etype);
  // \endcond
};
 
/** \brief Creates an empty graph for construction and loading
 */
Ngraph* createEmptyGraph();
 
/** \brief Cleans up the memory of the graph
 * \param g the graph
 */
void destroyGraph(Ngraph* g);

/** \brief Runs a validity check over the entire graph
 *
 */
bool checkValidity(Ngraph* g);

 void writeVTK(Ngraph* g,const char* prefix,GraphTag* tag=NULL,etype t=NO_TYPE);
 
} //namespace



#ifdef KOKKOS_ENABLED
#include <Kokkos_Core.hpp>
#define KOKKOS_FOR_VERTS(g,v)                                           \
  Kokkos::parallel_for(g->numLocalVtxs(),KOKKOS_LAMBDA(uintptr_t i) {   \
      agi::GraphVertex* v = reinterpret_cast<agi::GraphVertex*>((char*)(i+1));
      
#define KOKKOS_END_FOR() });
#endif


#endif
