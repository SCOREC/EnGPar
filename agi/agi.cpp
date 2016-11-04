#include "agi.h"
#include <cstdlib>
#include <stdint.h>
#include <iostream>
namespace agi {

Ngraph::Ngraph() {
  num_verts=0;
  num_types=0;
  weights=NULL;
  for (int i=0;i<MAX_TYPES;i++) {
    num_edges[i]=0;
    out_degree_list[i]=NULL;
    out_edges[i]=NULL;
  }
}

Ngraph::~Ngraph() {
  if (weights)
    delete [] weights;
  for (int i=0;i<MAX_TYPES;i++) {
    if (out_degree_list[i])
      delete [] out_degree_list[i];
    if (out_edges[i])
      delete [] out_edges[i];
  }
}

double Ngraph::weight(GraphVertex* vtx) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  return weights[index];
}
int  Ngraph::degree(GraphVertex* vtx,etype type) const {
  uintptr_t index =(uintptr_t)(vtx)-1;
  return out_degree_list[type][index+1]-out_degree_list[type][index];
}
  
GraphIterator* Ngraph::edges(GraphVertex* vtx,etype type) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  return reinterpret_cast<GraphIterator*>(out_edges[type]
                                          +out_degree_list[type][index]);
}

double Ngraph::weight(GraphEdge* edge) const {
  Edge* e = reinterpret_cast<Edge*>(edge);
  return e->weight;
}

GraphVertex* Ngraph::u(GraphEdge* edge) const {
  Edge* e = reinterpret_cast<Edge*>(edge);
  return reinterpret_cast<GraphVertex*>((char*)e->in+1);
}
GraphVertex* Ngraph::v(GraphEdge* edge) const {
  Edge* e = reinterpret_cast<Edge*>(edge);
  return reinterpret_cast<GraphVertex*>((char*)e->out+1);
}

GraphVertex* Ngraph::other(GraphEdge* edge,GraphVertex* vtx) const {
  if (isEqual(u(edge),vtx))
    return v(edge);
  return u(edge);
}


GraphIterator* Ngraph::begin(etype type) const {
  
  if (type == -1)
    return reinterpret_cast<GraphIterator*>((char*)1);
  return reinterpret_cast<GraphIterator*>(out_edges[type]);

}
#include <stdio.h>
GraphVertex* Ngraph::iterate(GraphIterator*& itr, GraphVertex*& vtx) const {
  uintptr_t index = (uintptr_t)(itr);
  if (index==num_verts+1)  {
    vtx=NULL;
    itr=NULL;
    return NULL;
  }
  vtx = reinterpret_cast<GraphVertex*>((int*)index);
  itr = reinterpret_cast<GraphIterator*>((char*)((uintptr_t)itr+1));
  return vtx;
}
GraphEdge* Ngraph::iterate(GraphIterator*& itr, GraphEdge*& edge) const {
  Edge* temp = reinterpret_cast<Edge*>(itr);
  edge = reinterpret_cast<GraphEdge*>(temp);
  temp++;
  itr=reinterpret_cast<GraphIterator*>(temp);
  return edge;
}

bool Ngraph::isEqual(GraphVertex* u,GraphVertex* v) const {
  return u==v;
}
  
//Protected functions

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
