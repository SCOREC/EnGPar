#ifndef AGI_H
#define AGI_H
#include <cstdlib>
#include <vector>
#include <map>
#include <cassert>
#include "Edge.h"

namespace agi {

//Definitions for edge types
#define MAX_TYPES 10 //static size of edge types
typedef int etype; 
#define SPLIT_TYPE 9 //predefined edge type for split vtx

  
class GraphVertex;
class GraphEdge;
class GraphIterator;

class Ngraph {

  //TODO: Make distributed version
  //TODO: Make ghost layer?

protected:
  // number of vertices and edges
  size_t num_verts;
  int num_types;
  size_t num_edges[MAX_TYPES];
  
  // vertex weights
  // size = num_verts
  double* weights;

  //TODO: add coordinate container
  
  /*
    offset list from a vertex to the out_vertices array
    degree of a vertex can be found by out_degree_list[i+1] - out_degree_list[i]
    beginning of edges can be found at out_vertices[out_degree_list[i]]
    size = num_verts+1
  */
  int* out_degree_list[MAX_TYPES];
  
  Edge* out_edges[MAX_TYPES];

  //TODO: do we want to have in edges as well?
  
public:
  Ngraph();
  virtual ~Ngraph();
  
  //Part Information
  size_t numVtxs() const {return num_verts;}
  int numEdgeTypes() const {return num_types;}
  size_t numEdges(etype i=0) const {return num_edges[i];}

  //Vertex Operations
  double weight(GraphVertex*) const;

  int owner(GraphVertex*) const {return 0;}

  const std::vector<double>& coordinates(GraphVertex*) const {};

  int degree(GraphVertex*,etype) const;
  
  GraphIterator* edges(GraphVertex*,etype) const;

  //Edge Operations
  double weight(GraphEdge*) const;

  GraphVertex* u(GraphEdge*) const;
  GraphVertex* v(GraphEdge*) const;

  GraphVertex* other(GraphEdge*,GraphVertex*) const;

  //Traversal
  GraphIterator* begin(etype i = -1) const;
  GraphVertex* iterate(GraphIterator*&,GraphVertex*&) const;
  GraphEdge* iterate(GraphIterator*&,GraphEdge*&) const;
  //virtual void destroy(GraphIterator*) = 0;

  //Utility
  bool isEqual(GraphVertex*,GraphVertex*) const;
  virtual void migrate(std::map<GraphVertex*,int>&) = 0;
  
 protected:
  etype addEdgeType() {assert(num_types!=MAX_TYPES);return num_types++;}
  Edge makeEdge(int in,int out, double w) {Edge e = {in,out,w}; return e;}
  
  /*
  void create_csr(int num_verts, int num_edges, int* srcs,
                  int* dsts, int* wgts);
  */
};

} //namespace

#endif
