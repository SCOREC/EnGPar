#include "ngraph.h"
#include "HyperEdgeIterator.h"
#include <cstdlib>
#include <stdint.h>
#include <iostream>
#include <PCU.h>
#include <cstring>
namespace agi {

Ngraph::Ngraph() {
  isHyperGraph=false;
  num_global_verts=0;
  num_local_verts=0;
  num_ghost_verts=0;
  
  num_types=0;
  local_weights=NULL;
  local_coords=NULL;
  for (int i=0;i<MAX_TYPES;i++) {
    edge_unmap[i] = NULL;
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

void Ngraph::constructGraph(bool isHG,
			    std::vector<gid_t>& verts,
			    std::vector<gid_t>& edge_ids,
			    std::vector<lid_t>& degs,
			    std::vector<gid_t>& pins_to_verts,
			    std::unordered_map<gid_t,part_t>& owns) {
  destroyData();
  isHyperGraph=isHG;
  num_local_verts=verts.size();
  local_unmap = new gid_t[num_local_verts];
  local_weights = NULL;
  local_coords = NULL;
  num_types = 0;
  etype t = addEdgeType();
  degree_list[t] = new lid_t[num_local_verts+1];
  degree_list[t][0] = 0;
  for (lid_t i=0;i<verts.size();i++) {
    local_unmap[i] = verts[i];
    vtx_mapping[verts[i]]=i;
    degree_list[t][i+1]=0;//degree_list[t][i]+degs[i];
    //set local_weights
    //set local coords
  }

  num_ghost_verts=0;
  num_local_edges[t] = edge_ids.size();
  num_local_pins[t] = pins_to_verts.size();
  edge_weights[t] = NULL;
  edge_unmap[t] = new gid_t[edge_ids.size()];
  pin_degree_list[t] = new lid_t[degs.size()+1];
  pin_degree_list[t][0]=0;
  pin_list[t] = new lid_t[pins_to_verts.size()];
  for (lid_t i=0;i<edge_ids.size();i++) {
    //set edge_weight
    gid_t gid = edge_ids[i];
    edge_mapping[t][gid]=i;
    edge_unmap[t][i]=gid;
    pin_degree_list[t][i+1]=pin_degree_list[t][i]+degs[i];
    for (lid_t j=pin_degree_list[t][i];j<pin_degree_list[t][i+1];j++) {
      gid_t v = pins_to_verts[j];
      map_t::iterator vitr = vtx_mapping.find(v);
      if (vitr!=vtx_mapping.end()&&vitr->second<num_local_verts) {
	if (!isHyperGraph&&!(j%2))
	  degree_list[t][vitr->second+1]++;
	else if (isHyperGraph)
	  degree_list[t][vitr->second+1]++;
	pin_list[t][j]=vitr->second;
      }
      else {
	if (vitr==vtx_mapping.end()) {
	  vtx_mapping[v]=num_local_verts+num_ghost_verts;
	  pin_list[t][j] = num_local_verts+num_ghost_verts++;
	}
	else {
	  pin_list[t][j]=vitr->second;
	}
      }
    }
  }
  if (!isHyperGraph)
    num_local_edges[t] = num_local_pins[t]/2;
  for (lid_t i=1;i<num_local_verts+1;i++) {
    degree_list[t][i]+=degree_list[t][i-1];
  }
  uint64_t* temp_counts = (uint64_t*)malloc(num_local_verts*sizeof(uint64_t));
  
  std::memcpy(temp_counts, degree_list[t], num_local_verts*sizeof(uint64_t));
  edge_list[t] = new lid_t[degree_list[t][num_local_verts]];
  ghost_unmap = new lid_t[num_ghost_verts];
  owners = new part_t[num_ghost_verts];
  for (lid_t i=0;i<edge_ids.size();i++) {
    for (lid_t j=pin_degree_list[t][i];j<pin_degree_list[t][i+1];j++) {
      lid_t u = pin_list[t][j];
      if (u>=num_local_verts){
	ghost_unmap[u-num_local_verts] = pins_to_verts[j];
	owners[u-num_local_verts] = owns[pins_to_verts[j]];
      }
      if (isHyperGraph) {
	if (u<num_local_verts)
	  edge_list[t][temp_counts[u]++] = i;
      }
      else {
	lid_t v = pin_list[t][++j];
	if (v>=num_local_verts){
	  ghost_unmap[v-num_local_verts] = pins_to_verts[j];
	  owners[v-num_local_verts] = owns[pins_to_verts[j]];
	}
	if (u<num_local_verts)
	  edge_list[t][temp_counts[u]++] = v;
      }
    }
  }
  free(temp_counts);
  //Setup global counters
  num_global_verts = PCU_Add_Long(num_local_verts);
  num_global_edges[t] = PCU_Add_Long(num_local_edges[t]);
  num_global_pins[t] = PCU_Add_Long(num_local_pins[t]);
}
Ngraph::~Ngraph() {
  destroyData();
}

void Ngraph::destroyData() {
  if (local_weights)
    delete [] local_weights;
  if (local_coords)
    delete [] local_coords;
  for (int i=0;i<MAX_TYPES;i++) {
    if (edge_unmap[i])
      delete [] edge_unmap[i];
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
  return edge_weights[type][edge_list[type][id]];
}

lid_t Ngraph::localID(GraphEdge* edge) const {
  uintptr_t id = (uintptr_t)(edge)-1;
  return id/num_types;
}

gid_t Ngraph::globalID(GraphEdge* edge) const {
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
    return 0;
  if (edge==NULL)
    return 0;
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  lid_t index = id;
  return pin_degree_list[type][index+1]-pin_degree_list[type][index];
}

PinIterator* Ngraph::pins(GraphEdge* edge) const {
  if (!isHyperGraph)
    return NULL;
  uintptr_t id = (uintptr_t)(edge)-1;
  etype type = id%num_types;
  id/=num_types;
  lid_t index = id;//edge_list[type][id];
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
  if (isHyperGraph)
    return new HyperEdgeIterator(t,num_types,num_local_edges[t]);
  return new EdgeIterator(t,num_types,0,num_local_edges[t]);
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
  itr->iterate();
  if (isHyperGraph&&!itr->isHyper()) {
    index-=1;
    etype t = index%num_types;
    index/=num_types;
    return (GraphEdge*)(num_types*edge_list[t][index]+t+1);
  }
  return (GraphEdge*)(index);
}
GraphVertex* Ngraph::iterate(PinIterator*& itr) const {
  lid_t* e = reinterpret_cast<lid_t*>(itr);
  uintptr_t id = *e+1;
  GraphVertex* vtx = reinterpret_cast<GraphVertex*>((char*)id);
  itr = reinterpret_cast<PinIterator*>(++e);
  return vtx;
}
GraphVertex* Ngraph::iterate(GraphIterator*& itr) const {
  GraphVertex* vtx;
  if (itr->isH) {
    if (itr->count >= degree(itr->edge)) {
      GraphEdge* edge = iterate(itr->eitr);
      if (!edge)
	return NULL;
      itr->setEdge(edge,pins(edge));
    }
    vtx = iterate(itr->pitr);
    itr->count++;
  }
  else {
    itr->edge = iterate(itr->eitr);
    vtx = v(itr->edge);
  }
  return vtx;
}
  
GraphEdge* Ngraph::edge(GraphIterator* itr) const {
  return itr->edge;
}
  
void Ngraph::destroy(EdgeIterator* itr) const {
  delete itr;
}
void Ngraph::destroy(GraphIterator* itr) const {
  delete itr->eitr;
  delete itr;
}

bool Ngraph::isEqual(GraphVertex* u,GraphVertex* v) const {
  return u==v;
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

void destroyGraph(Ngraph* g) {delete g;}
}
