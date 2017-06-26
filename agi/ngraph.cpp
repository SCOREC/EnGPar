#include "ngraph.h"
#include <cstdlib>
#include <stdint.h>
#include <iostream>
#include <cstring>
#include <engpar_support.h>
#include "Iterators/HyperEdgeIterator.h"
#include "Iterators/PinIterator.h"
#include "Iterators/GraphIterator.h"

namespace agi {

 
const wgt_t& Ngraph::weight(GraphVertex* vtx) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  if (index>=numTotalVtxs()){
    printf("[ERROR] invalid vertex given to weight(vtx)\n");
    throw 1;
  }    
  else if (index>=num_local_verts) {
    printf("[ERROR] weights unknown for ghost vertices\n");
    throw 2;
  }
  return local_weights[index];
}
const coord_t& Ngraph::coord(GraphVertex* vtx) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  if (index>=numTotalVtxs()){
    printf("[ERROR] invalid vertex given to coord(vtx)\n");
    throw 1;
  }    
  else if (index>=num_local_verts) {
    printf("[ERROR] coordinates unknown for ghost vertices\n");
    throw 2;
  }
  return local_coords[index];
}


int Ngraph::owner(GraphVertex* vtx) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  if (index>=num_local_verts+num_ghost_verts) {
    fprintf(stderr,"[ERROR] invalid vertex given to owner(vtx)\n");
    return -1;
  }
  if (index<num_local_verts)
    return PCU_Comm_Self();
  assert(PCU_Comm_Peers()>1);
  index-=num_local_verts;
  return owners[index];
}

part_t Ngraph::originalOwner(GraphVertex* vtx) const {
  assert(original_owners);
  uintptr_t index = (uintptr_t)(vtx)-1;
  if (index>=num_local_verts+num_ghost_verts) {
    fprintf(stderr,"[ERROR] invalid vertex given to owner(vtx)\n");
    return -1;
  }
  return original_owners[index];
}

void Ngraph::setOriginalOwners() {
  if (EnGPar_Is_Log_Open()) {
    char message[25];
    sprintf(message,"setOriginalOwners\n");
    EnGPar_Log_Function(message);
  }

  assert(!original_owners);
  original_owners = new part_t[num_local_verts];
  for (lid_t i=0;i<num_local_verts;i++)
    original_owners[i]=PCU_Comm_Self();
  if (EnGPar_Is_Log_Open()) {
    EnGPar_End_Function();
  }
}

void Ngraph::setOriginalOwners(std::vector<part_t>& oos) {
  assert(!original_owners);
  original_owners = new part_t[num_local_verts];
  for (lid_t i=0;i<num_local_verts;i++) 
    original_owners[i]=oos[i];
}

void Ngraph::getResidence(GraphEdge* e, Peers& residence) const {
  agi::PinIterator* pitr = pins(e);
  agi::GraphVertex* vtx;
  lid_t deg = degree(e);
  for (lid_t i=0;i<deg;i++) {
    vtx = iterate(pitr);
    residence.insert(owner(vtx));
  }
  destroy(pitr);
}
  
lid_t Ngraph::localID(GraphVertex* vtx) const {
  return  (uintptr_t)(vtx)-1;

}
gid_t Ngraph::globalID(GraphVertex* vtx) const {
  agi::lid_t lid= (uintptr_t)(vtx)-1;
  if (lid>=num_local_verts)
    return ghost_unmap[lid-num_local_verts];
  return  local_unmap[lid];
}

GraphVertex* Ngraph::find(GraphVertex* vtx) const {
  return findGID(globalID(vtx));
}

wgt_t Ngraph::weight(GraphEdge* edge) const {
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  if (isHyperGraph)
    return edge_weights[type][edge_list[type][id]];
  return edge_weights[type][id];
}

lid_t Ngraph::localID(GraphEdge* edge) const {
  uintptr_t id = (uintptr_t)(edge)-1;
  return id/num_types;
}

gid_t Ngraph::globalID(GraphEdge* edge) const {
  if (!isHyperGraph)
    return localID(edge);
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  return edge_unmap[type][id/num_types];
}

lid_t Ngraph::u(lid_t e, etype t) const {
  bool found = false;
  lid_t index = 0;
  lid_t bound_low=0;
  lid_t bound_high = numLocalVtxs();
  while (!found) {
    index = (bound_high+bound_low)/2;
    if (degree_list[t][index]<= e&& degree_list[t][index+1]>e) 
      found=true;
    else if (degree_list[t][index]<=e)
      bound_low=index;
    else
      bound_high=index;

  }
  return index;
}
GraphVertex* Ngraph::u(GraphEdge* edge) const {
  lid_t lid = (uintptr_t)(edge)-1;
  etype type = lid%num_types;
  lid/=num_types;
  lid_t vid = u(lid,type);
  return reinterpret_cast<GraphVertex*>(vid+1);
}

GraphVertex* Ngraph::v(GraphEdge* edge) const {
  if (isHyperGraph) {
    fprintf(stderr,"v(edge) not supported in hypergraph mode");
    return NULL;
  }
  if (edge==NULL) 
    return NULL;
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

GraphIterator* Ngraph::adjacent(GraphVertex* vtx, etype type) const {
  uintptr_t index = (uintptr_t)(vtx)-1;
  EdgeIterator* eitr = new EdgeIterator(type,num_types,(lid_t*)degree_list[type][index],degree(vtx,type));
  return new GraphIterator(eitr,isHyperGraph);
}


lid_t  Ngraph::degree(GraphEdge* edge) const {
  if (!isHyperGraph)
    return 2;
  if (edge==NULL)
    return 0;
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  lid_t index = id;
  return pin_degree_list[type][index+1]-pin_degree_list[type][index];
}

PinIterator* Ngraph::pins(GraphEdge* edge) const {
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  if (!isHyperGraph) {
    return new PinIterator(reinterpret_cast<lid_t*>(u(edge)),
                           reinterpret_cast<lid_t*>((char*)(edge_list[type][id]+1)));
  }
  return new PinIterator((pin_list[type]+pin_degree_list[type][id]),
                         pin_list[type]+pin_degree_list[type][id+1]);
}



bool Ngraph::isEqual(GraphVertex* u,GraphVertex* v) const {
  return u==v;
}

PartitionMap* Ngraph::getPartition() {
  if (EnGPar_Is_Log_Open()) {
    char message[20];
    sprintf(message,"getPartition\n");
    EnGPar_Log_Function(message);
  }

  PCU_Comm_Begin();
  VertexIterator* itr = begin();
  GraphVertex* vtx;
  while ((vtx = iterate(itr))) {
    gid_t v = globalID(vtx);
    PCU_COMM_PACK(originalOwner(vtx),v);
  }
  PCU_Comm_Send();
  PartitionMap* map = new PartitionMap;
  while (PCU_Comm_Receive()) {
    gid_t v;
    PCU_COMM_UNPACK(v);
    map->insert(std::make_pair(v,PCU_Comm_Sender()));
  }
  if (EnGPar_Is_Log_Open()) {
    EnGPar_End_Function();
  }
  PCU_Barrier();
  return map;
}

void Ngraph::setEdgeWeights(std::vector<wgt_t>& wgts, etype t) {
  assert(!edge_weights[t]);
  if (wgts.size()==0) {
    edge_weights[t] = new wgt_t[num_local_edges[t]];
    for (gid_t i=0;i<num_local_edges[t];i++)
      edge_weights[t][i]=1;
    return;
  }
  assert(wgts.size()==num_local_edges[t]);
  edge_weights[t] = new wgt_t[num_local_edges[t]];
  memcpy(edge_weights[t],&(wgts[0]),num_local_edges[t]*sizeof(wgt_t));
}

//Protected functions

void Ngraph::makeEdgeArray(etype t, int count) {
  edge_unmap[t] = new gid_t[count];
  edge_weights[t] = new wgt_t[count];
}

void Ngraph::setEdge(lid_t lid,gid_t gid, wgt_t w,etype t) {
  edge_unmap[t][lid] = gid;
  edge_weights[t][lid] = w;
  edge_mapping[t][gid]=lid;
  
}

void destroyGraph(Ngraph* g) {
  if (EnGPar_Is_Log_Open()) {
    char message[25];
    sprintf(message,"destroyGraph\n");
    EnGPar_Log_Function(message);
  }

  delete g;
  if (EnGPar_Is_Log_Open()) {
    EnGPar_End_Function();
  }  

}


}
