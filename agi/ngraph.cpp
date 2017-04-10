#include "ngraph.h"
#include "HyperEdgeIterator.h"
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

Ngraph::~Ngraph() {
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

void Ngraph::sendVertex(GraphVertex* vtx, part_t toSend) {
  //char vertex[100];
  gid_t gid = globalID(vtx);
  PCU_COMM_PACK(toSend,gid);
  GraphIterator* gitr = adjacent(vtx);
  GraphVertex* other;
  lid_t deg = degree(vtx);
  PCU_COMM_PACK(toSend,deg);
  GraphEdge* prevEdge=NULL;
  //sprintf(vertex,"gid: %lu\nHas degree: %lu\nEdges:\n",gid,deg); 
  while ((other = iterate(gitr))) {
    if (isHyperGraph) {
      GraphEdge* e = edge(gitr);
      if (prevEdge!=e) {
	prevEdge=e;
	gid_t edge_gid = globalID(e);
	PCU_COMM_PACK(toSend,edge_gid);
	//sprintf(vertex,"%s%lu\n",vertex,edge_gid);
      }
    }
    else {
      gid_t other_gid = globalID(other);
      PCU_COMM_PACK(toSend,other_gid);
      //sprintf(vertex,"%s%lu\n",vertex,other_gid);
    }
  }
  //printf("%s",vertex);
}

void Ngraph::recvVertex() {
  //char vertex[100];
  gid_t gid;
  PCU_COMM_UNPACK(gid);
  lid_t deg;
  PCU_COMM_UNPACK(deg);
  //sprintf(vertex,"gid: %lu\nHas degree: %lu\nEdges:\n",gid,deg);
  for (lid_t i=0; i<deg;i++) {
    if (isHyperGraph) {
      gid_t edge_gid;
      PCU_COMM_UNPACK(edge_gid);
      //sprintf(vertex,"%s%lu\n",vertex,edge_gid);
    }
    else {
      gid_t other_gid;
      PCU_COMM_UNPACK(other_gid);
      //sprintf(vertex,"%s%lu\n",vertex,other_gid);
    }
  }
  //printf("%s",vertex);
}

  
void Ngraph::migrate(Migration* plan) {
  Migration::iterator itr;
  PCU_Comm_Begin();
  for (itr = plan->begin();itr!=plan->end();itr++) {
    //printf("%d sending %p to %d\n",PCU_Comm_Self(),itr->first,itr->second);
    sendVertex(itr->first,itr->second);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    //printf("%d getting vertex\n",PCU_Comm_Self());
    recvVertex();
  }
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
