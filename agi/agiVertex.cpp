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

bool Ngraph::isEqual(GraphVertex* u,GraphVertex* v) const {
  return u==v;
}

}
