#include "ngraph.h"
#include <cstdlib>
#include <stdint.h>
#include <iostream>
namespace agi {

Ngraph::Ngraph() {
  num_local_verts=0;
  num_types=0;
  local_weights=NULL;
  for (int i=0;i<MAX_TYPES;i++) {
    num_local_edges[i]=0;
    degree_list[i]=NULL;
    edge_list[i]=NULL;
    pin_degree_list[i]=NULL;
    pin_list[i]=NULL;
  }
}

Ngraph::~Ngraph() {
  if (local_weights)
    delete [] local_weights;
  for (int i=0;i<MAX_TYPES;i++) {
    if (degree_list[i])
      delete [] degree_list[i];
    if (edge_list[i])
      delete [] edge_list[i];
    if (pin_degree_list[i])
      delete [] pin_degree_list[i];
    if (pin_list[i])
      delete [] pin_list[i];
  }
}

double Ngraph::weight(GraphVertex* vtx) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  return local_weights[index];
}

double Ngraph::weight(GraphEdge* edge) const {
  Edge* e = reinterpret_cast<Edge*>(edge);
  return e->weight;
}

lid_t  Ngraph::degree(GraphVertex* vtx,etype type) const {
  uintptr_t index =(uintptr_t)(vtx)-1;
  return degree_list[type][index+1]-degree_list[type][index];
}
  
EdgeIterator* Ngraph::edges(GraphVertex* vtx,etype type) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  EdgeIterator* eitr = new EdgeIterator(type,edge_list[type]
                                        + degree_list[type][index]);
  return eitr;
}

lid_t  Ngraph::degree(GraphEdge* edge) const {
  Edge* e = reinterpret_cast<Edge*>(edge);
  uintptr_t index = e->lid;
  etype type = e->type;
  return pin_degree_list[type][index+1]-pin_degree_list[type][index];
}

PinIterator* Ngraph::pins(GraphEdge* edge) const {
  Edge* e = reinterpret_cast<Edge*>(edge);
  lid_t index = e->lid;
  etype type = e->type;
  return reinterpret_cast<PinIterator*>(pin_list[type]+
                                        pin_degree_list[type][index]);
}


VertexIterator* Ngraph::begin() const {
  return reinterpret_cast<VertexIterator*>((char*)1);
}
  
GraphVertex* Ngraph::iterate(VertexIterator*& itr) const {
  uintptr_t index = (uintptr_t)(itr);
  if (index==num_local_verts+1)  {
    itr=NULL;
    return NULL;
  }
  GraphVertex* vtx = reinterpret_cast<GraphVertex*>((lid_t*)index);
  itr = reinterpret_cast<VertexIterator*>((char*)(index+1));
  return vtx;
}
GraphEdge* Ngraph::iterate(EdgeIterator* itr) const {
  return reinterpret_cast<GraphEdge*>(es[itr->type]+*(itr->loc++));
}
GraphVertex* Ngraph::iterate(PinIterator*& itr) const {
  lid_t* e = reinterpret_cast<lid_t*>(itr);
  GraphVertex* vtx = reinterpret_cast<GraphVertex*>((char*)((*e)+1));
  itr = reinterpret_cast<PinIterator*>(++e);
  return vtx;
}

void Ngraph::destroy(EdgeIterator* itr) const {
  delete itr;
}

bool Ngraph::isEqual(GraphVertex* u,GraphVertex* v) const {
  return u==v;
}
  
//Protected functions

void Ngraph::makeEdgeArray(etype t, int count) {
  es[t] = new Edge[count];
}

void Ngraph::setEdge(lid_t lid,gid_t gid, wgt_t w,etype t) {
  es[t][lid] = Edge(lid,gid,w,t);
  edge_mapping[t][gid]=lid;
  
}
  /*
void Ngraph::create_csr(int nv, int ne, int* srcs,
                          int* dsts, int* wgts) {
  num_verts = nv;
  num_edges = ne;
  weights = new double[num_verts];
  out_vertices = new int[num_edges];
  out_weights = new double[num_edges];
  out_degree_list = new int[num_verts+1];

  for (size_t i = 0; i < num_edges; ++i)
    out_vertices[i] = 0;
  for (size_t i = 0; i < num_edges; ++i)
    out_weights[i] = 0;
  for (size_t i = 0; i < num_verts+1; ++i)
    out_degree_list[i] = 0;

  int* temp_counts = new int[num_verts];
  for (size_t i = 0; i < num_verts; ++i)
    temp_counts[i] = 0;
  for (size_t i = 0; i < num_edges; ++i)
    ++temp_counts[srcs[i]];
  for (size_t i = 0; i < num_verts; ++i)
    out_degree_list[i+1] = out_degree_list[i] + temp_counts[i];
  std::copy(out_degree_list, out_degree_list + num_verts, temp_counts);
  for (size_t i = 0; i < num_edges; ++i) {
    out_vertices[temp_counts[srcs[i]]] = dsts[i];
    out_weights[temp_counts[srcs[i]]++] = wgts[i];
  }
  
  delete [] temp_counts;

}
  */
}
