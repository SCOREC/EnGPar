#ifndef NGRAPH_H__
#define NGRAPH_H__
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <map>
#include <cassert>
#include "Edge.h"
#include "agi.h"
#include "EdgeIterator.h"

namespace agi {

class GraphVertex;
class GraphEdge;
class VertexIterator;
class PinIterator;
class Ngraph {
  //TODO: Try to compress global id
protected:

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

  //edge weights
  // size=num_edges
  //TODO: think of a way to use struct of arrays instead
  //TODO: allow disable hypergraph for simple graphs to minimize memory
  gid_t* edge_ids[MAX_TYPES];
  wgt_t* edge_weights[MAX_TYPES];
  
  //Edge* es[MAX_TYPES];
  
  //TODO: do we want ghost weights?
  //wgt_t* ghost_weights

  //TODO: add coordinate container
  //coord_t* local_coords;
  
  lid_t* degree_list[MAX_TYPES];
  lid_t* edge_list[MAX_TYPES];
  lid_t* pin_degree_list[MAX_TYPES];
  lid_t* pin_list[MAX_TYPES];

  //TODO: discuss using C++11 to get unordered map
  typedef std::unordered_map<gid_t,lid_t> map_t;
  map_t vtx_mapping;
  map_t edge_mapping[MAX_TYPES];
  //TODO: Tack ghost unmap on top of local_unmap
  gid_t* local_unmap;
  gid_t* ghost_unmap;
  part_t* owners;
  
public:
  Ngraph();
  ~Ngraph();

  //Global Part Information
  gid_t numGlobalVtxs() const {return num_global_verts;}
  gid_t numGlobalEdges(etype i=0) const {return num_global_edges[i];}
  gid_t numGlobalPins(etype i=0) const {return num_global_pins[i];}
  
  //Local Part Information
  lid_t numLocalVtxs() const {return num_local_verts;}
  lid_t numGhostVtxs() const {return num_ghost_verts;}
  lid_t numTotalVtxs() const {return num_local_verts+num_ghost_verts;}
  int numEdgeTypes() const {return num_types;}
  lid_t numLocalEdges(etype i=0) const {return num_local_edges[i];}
  lid_t numLocalPins(etype i=0) const {return num_local_pins[i];}
  
  //Vertex Operations
  wgt_t weight(GraphVertex*) const;
  part_t owner(GraphVertex*) const;
  const std::vector<double>& coordinates(GraphVertex*) const {};
  lid_t localID(GraphVertex*) const;
  gid_t globalID(GraphVertex*) const;
  GraphVertex* find(GraphVertex* vtx) const;

  //Edge Operations
  double weight(GraphEdge*) const;
  lid_t u(lid_t) const;
  GraphVertex* v(GraphEdge*) const;
  
  //Adjacency Operations
  lid_t degree(GraphVertex*,etype t=0) const;
  EdgeIterator* edges(GraphVertex*,etype t=0) const;
  lid_t degree(GraphEdge*) const;
  PinIterator* pins(GraphEdge*) const;
  
  //Iterator Traversal
  VertexIterator* begin() const;
  GraphVertex* findGID(gid_t gid) const;
  EdgeIterator* begin(etype t) const;
  //Iterate through vertices
  GraphVertex* iterate(VertexIterator*&) const;
  //Iterate through Edges
  GraphEdge* iterate(EdgeIterator*&) const;
  //Iterate through Pins
  GraphVertex* iterate(PinIterator*&) const;
  //Destroys iterator
  void destroy(EdgeIterator*) const;
  //Utility
  bool isEqual(GraphVertex*,GraphVertex*) const;
  virtual void migrate(std::map<GraphVertex*,int>&) = 0;
  
 protected:
  etype addEdgeType() {assert(num_types!=MAX_TYPES);return num_types++;}

  void makeEdgeArray(etype t,int count);
  void setEdge(lid_t,gid_t,wgt_t,etype);
  
  /*
  void create_csr(int num_verts, int num_edges, int* srcs,
                  int* dsts, int* wgts);
  */
};

} //namespace

#endif
