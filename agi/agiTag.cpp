#include "ngraph.h"
#include "agiInfo.h"
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
  int Ngraph::getIntTag(GraphTag* t,GraphVertex* vtx) {
    uintptr_t index = (uintptr_t)(vtx)-1;
    return reinterpret_cast<Info*>(t)->getInt(index);
  }
  int Ngraph::getIntTag(GraphTag* t,GraphEdge* edge) {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    return reinterpret_cast<Info*>(t)->getInt(id);
  }
  double Ngraph::getDoubleTag(GraphTag* t,GraphVertex* vtx) {
    uintptr_t index = (uintptr_t)(vtx)-1;
    return reinterpret_cast<Info*>(t)->getDouble(index);
  }
  double Ngraph::getDoubleTag(GraphTag* t,GraphEdge* edge) {
    uintptr_t id = (uintptr_t)(edge)-1;
    id/=num_types;
    return reinterpret_cast<Info*>(t)->getDouble(id);
  }
  long Ngraph::getLongTag(GraphTag* t,GraphVertex* vtx) {
    uintptr_t index = (uintptr_t)(vtx)-1;
    return reinterpret_cast<Info*>(t)->getLong(index);
  }
  long Ngraph::getLongTag(GraphTag* t,GraphEdge* edge) {
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

}
