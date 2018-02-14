#include "ngraph.h"
#include "Iterators/HyperEdgeIterator.h"
#include "Iterators/PinIterator.h"
#include "Iterators/GraphIterator.h"

namespace agi {
  VertexIterator* Ngraph::begin() const {
    return reinterpret_cast<VertexIterator*>((char*)1);
  }
  GhostIterator* Ngraph::beginGhosts() const {
    return reinterpret_cast<GhostIterator*>((char*)num_local_verts+1);
  }
  GraphVertex* Ngraph::findGID(gid_t gid) const {
    return reinterpret_cast<GraphVertex*>((char*)(vtx_mapping.find(gid)->second)+1);
  }

  EdgeIterator* Ngraph::begin(etype t) const {
    if (isHyperGraph)
      return new HyperEdgeIterator(t,num_types,num_local_edges[t]);
    return new EdgeIterator(t,num_types,0,num_local_edges[t]);
  }
  
  GraphVertex* Ngraph::iterate(VertexIterator*& itr) const {
    lid_t index = (uintptr_t)(itr);
    if (index==num_local_verts+1)  {
      itr=NULL;
      return NULL;
    }
    GraphVertex* vtx = reinterpret_cast<GraphVertex*>((lid_t*)index);
    itr = reinterpret_cast<VertexIterator*>((char*)(index+1));
    return vtx;
  }
  GraphVertex* Ngraph::iterate(GhostIterator*& itr) const {
    lid_t index = (uintptr_t)(itr);
    if (index==num_local_verts+num_ghost_verts+1)  {
      itr=NULL;
      return NULL;
    }
    GraphVertex* vtx = reinterpret_cast<GraphVertex*>((lid_t*)index);
    itr = reinterpret_cast<GhostIterator*>((char*)(index+1));
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
    if (isHyperGraph) {
      if (itr->loc==itr->end)
        return NULL;
      lid_t* e = itr->loc;
      uintptr_t id = *e+1;
      GraphVertex* vtx = reinterpret_cast<GraphVertex*>((char*)id);
      itr->iterate();
      return vtx;
    }
    else {
      GraphVertex* vtx = reinterpret_cast<GraphVertex*>(itr->loc);
      if (itr->loc==itr->end)
        itr->end=NULL;
      itr->loc=itr->end;
      return vtx;
    }
  }
  GraphVertex* Ngraph::iterate(GraphIterator*& itr) const {
    GraphVertex* vtx;
    if (itr->isH) {
      if (itr->count >= degree(itr->edge)) {
        destroy(itr->pitr);
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
  void Ngraph::destroy(PinIterator* itr) const {
    delete itr;
  }
  void Ngraph::destroy(GraphIterator* itr) const {
    delete itr->eitr;
    delete itr;
  }
}
