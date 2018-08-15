#include "ngraph.h"
#include "agiInfo.h"
#include <PCU.h>
namespace agi {


  GraphTag* Ngraph::createIntTag(etype t) {
    lid_t size;
    int bytes = sizeof(int);
    if (t==VTX_TYPE)
      size = num_local_verts;
    else
      size = num_local_edges[t];
    Info* i = new Info(size,bytes);
    return reinterpret_cast<GraphTag*>(i);
  }
  GraphTag* Ngraph::createDoubleTag(etype t) {
    lid_t size;
    int bytes = sizeof(double);
    if (t==VTX_TYPE)
      size = num_local_verts;
    else
      size = num_local_edges[t];
    Info* i = new Info(size,bytes);
    return reinterpret_cast<GraphTag*>(i);
  }
  GraphTag* Ngraph::createLongTag(etype t) {
    lid_t size;
    int bytes = sizeof(long);
    if (t==VTX_TYPE)
      size = num_local_verts;
    else
      size = num_local_edges[t];
    Info* i = new Info(size,bytes);
    return reinterpret_cast<GraphTag*>(i);
  }
  void Ngraph::destroyTag(GraphTag* t) {
    delete reinterpret_cast<Info*>(t);
  }
  int Ngraph::getIntTag(GraphTag* t,GraphVertex* vtx) const {
    uintptr_t index = (uintptr_t)(vtx)-1;
    if (index>=(unsigned)num_local_verts)
      index-=num_local_verts;
    return reinterpret_cast<Info*>(t)->getInt(index);
  }
  int Ngraph::getIntTag(GraphTag* t,GraphEdge* edge) const {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    return reinterpret_cast<Info*>(t)->getInt(id);
  }
  double Ngraph::getDoubleTag(GraphTag* t,GraphVertex* vtx) const {
    uintptr_t index = (uintptr_t)(vtx)-1;
    if (index>=(unsigned)num_local_verts)
      index-=num_local_verts;
    return reinterpret_cast<Info*>(t)->getDouble(index);
  }
  double Ngraph::getDoubleTag(GraphTag* t,GraphEdge* edge) const {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    return reinterpret_cast<Info*>(t)->getDouble(id);
  }
  long Ngraph::getLongTag(GraphTag* t,GraphVertex* vtx) const {
    uintptr_t index = (uintptr_t)(vtx)-1;
    if (index >= (unsigned)num_local_verts)
      index -= num_local_verts;
    return reinterpret_cast<Info*>(t)->getLong(index);
  }
  long Ngraph::getLongTag(GraphTag* t,GraphEdge* edge) const {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    return reinterpret_cast<Info*>(t)->getLong(id);
  }
  void Ngraph::setIntTag(GraphTag* t,GraphVertex* vtx,int val) {
    uintptr_t index = (uintptr_t)(vtx)-1;
    reinterpret_cast<Info*>(t)->setInt(index,val);
  }
  void Ngraph::setIntTag(GraphTag* t,GraphEdge* edge ,int val) {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    reinterpret_cast<Info*>(t)->setInt(id,val);
  }
  void Ngraph::setDoubleTag(GraphTag* t,GraphVertex* vtx,double val){
    uintptr_t index = (uintptr_t)(vtx)-1;
    reinterpret_cast<Info*>(t)->setDouble(index,val);
  }
  void Ngraph::setDoubleTag(GraphTag* t,GraphEdge* edge,double val) {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    reinterpret_cast<Info*>(t)->setDouble(id,val);
  }
  void Ngraph::setLongTag(GraphTag* t,GraphVertex* vtx,long val) {
    uintptr_t index = (uintptr_t)(vtx)-1;
    reinterpret_cast<Info*>(t)->setLong(index,val);
  }
  void Ngraph::setLongTag(GraphTag* t,GraphEdge* edge,long val) {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    reinterpret_cast<Info*>(t)->setLong(id,val);
  }
  GraphTag* Ngraph::createIntGhostTag(IntGhostTag callback) {
    lid_t size = num_ghost_verts;
    int bytes = sizeof(int);
    Info* i = new Info(size, bytes);
    PCU_Comm_Begin();
    GraphVertex* v;
    GhostIterator* vitr = beginGhosts();
    while ((v = iterate(vitr))) {
      gid_t gid = globalID(v);
      PCU_COMM_PACK(owner(v),gid);
    }
    PCU_Comm_Send();
    std::vector<std::pair<part_t, std::pair<gid_t, int> > >data;
    while (PCU_Comm_Receive()) {
      part_t sender = PCU_Comm_Sender();
      agi::gid_t gid;
      PCU_COMM_UNPACK(gid);
      v = getVertex(vtx_mapping[gid]);
      int val = callback(this,v);
      data.push_back(std::make_pair(sender, std::make_pair(gid, val)));
    }

    PCU_Comm_Begin();
    for (size_t i = 0; i < data.size(); i++) {
      PCU_COMM_PACK(data[i].first,data[i].second.first);
      PCU_COMM_PACK(data[i].first,data[i].second.second);
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      gid_t gid;
      int val;
      PCU_COMM_UNPACK(gid);
      PCU_COMM_UNPACK(val);
      v = getVertex(vtx_mapping[gid]);
      i->setInt(localID(v)-num_local_verts, val);
    }
    return reinterpret_cast<GraphTag*>(i);
  }
  GraphTag* Ngraph::createLongGhostTag(LongGhostTag callback) {
    lid_t size = num_ghost_verts;
    int bytes = sizeof(long);
    Info* i = new Info(size, bytes);
    PCU_Comm_Begin();
    GraphVertex* v;
    GhostIterator* vitr = beginGhosts();
    while ((v = iterate(vitr))) {
      gid_t gid = globalID(v);
      PCU_COMM_PACK(owner(v),gid);
    }
    PCU_Comm_Send();
    std::vector<std::pair<part_t, std::pair<gid_t, long> > >data;
    while (PCU_Comm_Receive()) {
      part_t sender = PCU_Comm_Sender();
      agi::gid_t gid;
      PCU_COMM_UNPACK(gid);
      v = getVertex(vtx_mapping[gid]);
      long val = callback(this,v);
      data.push_back(std::make_pair(sender, std::make_pair(gid, val)));
    }

    PCU_Comm_Begin();
    for (size_t i = 0; i < data.size(); i++) {
      PCU_COMM_PACK(data[i].first,data[i].second.first);
      PCU_COMM_PACK(data[i].first,data[i].second.second);
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      gid_t gid;
      long val;
      PCU_COMM_UNPACK(gid);
      PCU_COMM_UNPACK(val);
      v = getVertex(vtx_mapping[gid]);
      i->setLong(localID(v)-num_local_verts, val);
    }
    return reinterpret_cast<GraphTag*>(i);
  }

  GraphTag* Ngraph::createDoubleGhostTag(DoubleGhostTag callback) {
    lid_t size = num_ghost_verts;
    int bytes = sizeof(double);
    Info* i = new Info(size, bytes);
    PCU_Comm_Begin();
    GraphVertex* v;
    GhostIterator* vitr = beginGhosts();
    while ((v = iterate(vitr))) {
      gid_t gid = globalID(v);
      PCU_COMM_PACK(owner(v),gid);
    }
    PCU_Comm_Send();
    std::vector<std::pair<part_t, std::pair<gid_t, double> > >data;
    while (PCU_Comm_Receive()) {
      part_t sender = PCU_Comm_Sender();
      agi::gid_t gid;
      PCU_COMM_UNPACK(gid);
      v = getVertex(vtx_mapping[gid]);
      double val = callback(this,v);
      data.push_back(std::make_pair(sender, std::make_pair(gid, val)));
    }
    PCU_Comm_Begin();
    for (size_t i = 0; i < data.size(); i++) {
      PCU_COMM_PACK(data[i].first,data[i].second.first);
      PCU_COMM_PACK(data[i].first,data[i].second.second);
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      gid_t gid;
      double val;
      PCU_COMM_UNPACK(gid);
      PCU_COMM_UNPACK(val);
      v = getVertex(vtx_mapping[gid]);
      i->setDouble(localID(v)-num_local_verts, val);
    }
    return reinterpret_cast<GraphTag*>(i);

  }
}
