#include "ngraph.h"
#include <cstdlib>
#include <stdint.h>
#include <iostream>
#include <cstring>
#include <engpar_support.h>
#include "Iterators/HyperEdgeIterator.h"
#include "Iterators/PinIterator.h"
#include "Iterators/GraphIterator.h"
#include "agi_typeconvert.h"

namespace agi {

  GraphVertex* PNgraph::getVertex(lid_t lid) {
    return reinterpret_cast<GraphVertex*>( toPtr(lid+1) );
  }

const wgt_t& Ngraph::weight(GraphVertex* vtx) const {
  lid_t index = fromPtr(vtx);
  if (index>=numTotalVtxs()){
    EnGPar_Error_Message("invalid vertex given to weight(vtx)\n");
    throw 1;
  }    
  else if (index>=num_local_verts) {
    EnGPar_Error_Message("weights unknown for ghost vertices\n");
    throw 2;
  }
  return local_weights[index];
}
void Ngraph::setWeights(wgt_t* wgts) {
  if (!local_weights) {
    local_weights = new wgt_t[num_local_verts];
  }
  memcpy(local_weights,wgts,num_local_verts*sizeof(wgt_t));
}
bool Ngraph::hasCoords() const {
  return local_coords!=NULL;
}

void Ngraph::setCoords(coord_t* cs) {
  local_coords = new coord_t[num_local_verts];
  for (lid_t i=0;i<num_local_verts;i++)
    for (int j=0;j<3;j++)
      local_coords[i][j] = cs[i][j];
}
  
const coord_t& Ngraph::coord(GraphVertex* vtx) const {
  lid_t index = fromPtr(vtx);
  if (index>=numTotalVtxs()){
    EnGPar_Error_Message("invalid vertex given to coord(vtx)\n");
    throw 1;
  }    
  else if (index>=num_local_verts) {
    EnGPar_Error_Message("coordinates unknown for ghost vertices\n");
    throw 2;
  }
  return local_coords[index];
}


int Ngraph::owner(GraphVertex* vtx) const {
  lid_t index = fromPtr(vtx);
  if (index>=num_local_verts+num_ghost_verts) {
    EnGPar_Error_Message("invalid vertex given to owner(vtx)\n");
    return -1;
  }
  if (index<num_local_verts)
    return PCU_Comm_Self();
  index-=num_local_verts;
  return owners[index];
}

part_t Ngraph::originalOwner(GraphVertex* vtx) const {
  assert(original_owners);
  lid_t index = fromPtr(vtx);
  if (index>=num_local_verts+num_ghost_verts) {
    EnGPar_Error_Message("invalid vertex given to owner(vtx)\n");
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
  if (original_owners)
    return;
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

void Ngraph::resetOwnership() {
  if (EnGPar_Is_Log_Open()) {
    char message[25];
    sprintf(message,"resetOriginalOwners\n");
    EnGPar_Log_Function(message);
  }
  if (original_owners)
    delete [] original_owners;
  original_owners = NULL;
  setOriginalOwners();
  if (EnGPar_Is_Log_Open()) {
    EnGPar_End_Function();
  }
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
  lid_t index = fromPtr(vtx);
  return degree_list[type][index+1]-degree_list[type][index];
}
  
EdgeIterator* Ngraph::edges(GraphVertex* vtx,etype type) const {
  lid_t index = fromPtr(vtx);
  lid_t* foo = (lid_t*) toPtr(degree_list[type][index]);
  EdgeIterator* eitr = new EdgeIterator(type,num_types,foo,degree(vtx,type));
  return eitr;
}

GraphIterator* Ngraph::adjacent(GraphVertex* vtx, etype type) const {
  lid_t index = fromPtr(vtx);
  lid_t* foo = (lid_t*) toPtr(degree_list[type][index]);
  EdgeIterator* eitr = new EdgeIterator(type,num_types,foo,degree(vtx,type));
  return new GraphIterator(eitr,isHyperGraph);
}

bool Ngraph::isEqual(GraphVertex* u,GraphVertex* v) const {
  return u==v;
}

  void Ngraph::changeOwners(int* newRanks) {
    for (int i=0;i<num_ghost_verts;++i)
      owners[i] = newRanks[owners[i]];
  }

  VEVIterator* Ngraph::vev_begin(GraphVertex* vtx, etype t) const {
    lid_t index = fromPtr(vtx);
    lid_t start = vev_offsets[t][index]*num_types+t;
    return reinterpret_cast<VEVIterator*>(toPtr(start+1));
  }
  VEVIterator* Ngraph::vev_end(GraphVertex* vtx, etype t) const {
    lid_t index = fromPtr(vtx);
    lid_t start = vev_offsets[t][index+1]*num_types+t;
    return reinterpret_cast<VEVIterator*>(toPtr(start+1));
  }
  
}

