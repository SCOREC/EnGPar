
#ifndef AGI_H
#define AGI_H

#include <vector>
#include <map>

namespace agi {

class GraphVertex;
class GraphEdge;
class VertexIterator;
//TODO: do we need an EdgeIterator?

typedef int etype;

class Ngraph {

  //TODO: Make distributed version
  //TODO: Make ghost layer?

protected:
  //a flag controlling whether the graph is directed or undirected
  bool isUnd;
  
  // number of vertices and edges
  int num_verts;
  int num_edges;
  
  // vertex weights
  // size = num_verts
  int* weights;

  //TODO: add coordinate container
  
  /*
    offset list from a vertex to the out_vertices array
    degree of a vertex can be found by out_degree_list[i+1] - out_degree_list[i]
    beginning of edges can be found at out_vertices[out_degree_list[i]]
    size = num_verts+1
  */
  int* out_degree_list;

  // List of edges separated by the out_degree_list
  // size = num_edges
  int* out_vertices;

  //edge weights
  // size = num_edges
  int* out_weights;

  
public:
  Ngraph(bool isUndirected=true);
  virtual ~Ngraph();
  
  //Part Information
  
  int numVtx() const {return num_verts;}
  int numEdges() {return num_edges;}
    
  //Vertex Operations
  virtual double weight(GraphVertex*) const =0;

  virtual int owner(GraphVertex*) const =0;

  virtual const std::vector<double>& coordinates(GraphVertex*) const=0;

  virtual int degree(GraphVertex*,etype) const =0;
    
  virtual int* edges(GraphVertex*,etype) const =0;

  //Edge Operations
  virtual double weight(GraphEdge*) const =0;

  virtual GraphVertex* u(GraphEdge*) const =0;
  virtual GraphVertex* v(GraphEdge*) const =0;

  virtual GraphVertex* other(GraphEdge*,GraphVertex*) const =0;

  //Traversal
  virtual VertexIterator* begin() = 0;
  virtual VertexIterator* iterate(VertexIterator*) = 0;
  virtual void destroy(VertexIterator*) = 0;

  //Utility
  virtual bool isEqual(GraphVertex*,GraphVertex*) =0;
  virtual void migrate(std::map<GraphVertex*,int>&) = 0;
 protected:
  int out_deg(int n) const { return out_degree_list[n+1] - out_degree_list[n];}
  int* out_verts(int n) const {return out_vertices+out_degree_list[n];}
  int* out_wgts(int n) const {return out_weights+out_degree_list[n];}

  void create_csr(int num_verts, int num_edges, int* srcs,
                  int* dsts, int* wgts);
};

} //namespace

#endif
