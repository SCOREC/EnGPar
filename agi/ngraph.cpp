#include "ngraph.h"
#include <cstdlib>
#include <stdint.h>
#include <iostream>
#include <PCU.h>
namespace agi {

Ngraph::Ngraph() {
  isHyperGraph=false;
  num_global_verts=0;
  num_local_verts=0;
  num_ghost_verts=0;
  
  num_types=0;
  local_weights=NULL;

  for (int i=0;i<MAX_TYPES;i++) {
    edge_ids[i] = NULL;
    edge_weights[i] = NULL;
    num_global_edges[i]=0;
    num_local_edges[i]=0;
    degree_list[i]=NULL;
    edge_list[i]=NULL;
    num_global_pins[i]=0;
    num_local_pins[i]=0;
    pin_degree_list[i]=NULL;
    pin_list[i]=NULL;
  }

  local_unmap = NULL;
  ghost_unmap=NULL;
  owners = NULL;
}

Ngraph::~Ngraph() {
  if (local_weights)
    delete [] local_weights;
  
  for (int i=0;i<MAX_TYPES;i++) {
    if (edge_ids[i])
      delete [] edge_ids[i];
    if (edge_weights[i])
      delete [] edge_weights[i];
    if (degree_list[i])
      delete [] degree_list[i];
    if (edge_list[i])
      delete [] edge_list[i];
    if (pin_degree_list[i])
      delete [] pin_degree_list[i];
    if (pin_list[i])
      delete [] pin_list[i];
  }
  if (local_unmap)
    delete [] local_unmap;
  if (ghost_unmap)
    delete [] ghost_unmap;
  if (owners)
    delete [] owners;
}
 
double Ngraph::weight(GraphVertex* vtx) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  if (index>=num_local_verts) {
    printf("[ERROR] weights unknown for ghost vertices\n");
    return 0;
  }
  return local_weights[index];
}


int Ngraph::owner(GraphVertex* vtx) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  if (index>=num_local_verts+num_ghost_verts) {
    fprintf(stderr,"[ERROR] invalid vertex given to owner(vtx)\n");
    fprintf(stderr,"  Error vertex number: %d, Total vertices: %d\n",index,num_local_verts+num_ghost_verts);
    return -1;
  }
  if (index<num_local_verts)
    return PCU_Comm_Self();
  index-=num_local_verts;
  return owners[index];
}

lid_t Ngraph::localID(GraphVertex* vtx) const {
  return  (uintptr_t)(vtx)-1;

}
gid_t Ngraph::globalID(GraphVertex* vtx) const {
  return  local_unmap[(uintptr_t)(vtx)-1];
}

GraphVertex* Ngraph::find(GraphVertex* vtx) const {
  return findGID(globalID(vtx));
}


wgt_t Ngraph::weight(GraphEdge* edge) const {
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  return edge_weights[type][edge_list[type][id]];
}

lid_t Ngraph::u(lid_t e) const {
  bool found = false;
  lid_t index = 0;
  lid_t bound_low=0;
  lid_t bound_high = numLocalVtxs();
  while (!found) {
    index = (bound_high+bound_low)/2;
    if (degree_list[0][index]<= e&& degree_list[0][index+1]>e) 
      found=true;
    else if (degree_list[0][index]<=e)
      bound_low=index;
    else
      bound_high=index;

  }
  return index;
}

GraphVertex* Ngraph::v(GraphEdge* edge) const {
  if (isHyperGraph) {
    fprintf(stderr,"v(edge) not supported in hypergraph mode");
    return NULL;
  }
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  return reinterpret_cast<GraphVertex*>(edge_list[type][id]+1);
}


lid_t  Ngraph::degree(GraphVertex* vtx,etype type) const {
  uintptr_t index =(uintptr_t)(vtx)-1;
  return degree_list[type][index+1]-degree_list[type][index];
}
  
EdgeIterator* Ngraph::edges(GraphVertex* vtx,etype type) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  EdgeIterator* eitr = new EdgeIterator(type,num_types,(lid_t*)degree_list[type][index],degree(vtx,type));
  return eitr;
}

lid_t  Ngraph::degree(GraphEdge* edge) const {
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  lid_t index = edge_list[type][id];
  return pin_degree_list[type][index+1]-pin_degree_list[type][index];
}

PinIterator* Ngraph::pins(GraphEdge* edge) const {
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  lid_t index = edge_list[type][id];
  return reinterpret_cast<PinIterator*>(pin_list[type]+
                                        pin_degree_list[type][index]);
}


VertexIterator* Ngraph::begin() const {
  return reinterpret_cast<VertexIterator*>((char*)1);
}
GraphVertex* Ngraph::findGID(gid_t gid) const {
  return reinterpret_cast<GraphVertex*>((char*)(vtx_mapping.find(gid)->second));
}

EdgeIterator* Ngraph::begin(etype t) const {
  return new EdgeIterator(t,num_types,edge_list[t],num_local_edges[t]);
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
GraphEdge* Ngraph::iterate(EdgeIterator*& itr) const {
  if (itr->loc>=itr->end)
    return NULL;
  uintptr_t index = (uintptr_t)itr->loc;
  itr->loc = (gid_t*)(index+num_types);
  return (GraphEdge*)(index);
}
GraphVertex* Ngraph::iterate(PinIterator*& itr) const {
  lid_t* e = reinterpret_cast<lid_t*>(itr);
  uintptr_t id = *e+1;
  GraphVertex* vtx = reinterpret_cast<GraphVertex*>((char*)id);
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
  edge_ids[t] = new gid_t[count];
  edge_weights[t] = new wgt_t[count];
}

void Ngraph::setEdge(lid_t lid,gid_t gid, wgt_t w,etype t) {
  edge_ids[t][lid] = gid;
  edge_weights[t][lid] = w;
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
