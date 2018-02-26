#include "ngraph.h"
#include <stdint.h>
#include <iostream>
#include <cstring>
#include <engpar_support.h>
#include "Iterators/PinIterator.h"
#include "Iterators/GraphIterator.h"

namespace agi {

  void Ngraph::getResidence(GraphEdge* e, Peers& residence) const {
    //residence.clear();
    agi::PinIterator* pitr = pins(e);
    agi::GraphVertex* vtx;
    lid_t deg = degree(e);
    for (lid_t i=0;i<deg;i++) {
      vtx = iterate(pitr);
      residence.insert(owner(vtx));
    }
    destroy(pitr);
  }
  bool Ngraph::isResidentOn(GraphEdge* e,part_t peer) const {
    agi::PinIterator* pitr = pins(e);
    agi::GraphVertex* vtx;
    lid_t deg = degree(e);
    for (lid_t i=0;i<deg;i++) {
      vtx = iterate(pitr);
      if (owner(vtx)==peer) {
        destroy(pitr);
        return true;
      }
    }
    destroy(pitr);
    return false;
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

  lid_t  Ngraph::degree(GraphEdge* edge) const {
    if (!isHyperGraph)
      return 2;
    //This check is used for GraphIterator
    if (edge==NULL)
      return 0;
    uintptr_t id = (uintptr_t)(edge)-1;
    etype type = id%num_types;
    id/=num_types;
    return pin_degree_list[type][id+1]-pin_degree_list[type][id];
  }

  PinIterator* Ngraph::pins(GraphEdge* edge) const {
    uintptr_t id = (uintptr_t)(edge)-1;
    etype type = id%num_types;
    id/=num_types;
    if (!isHyperGraph) {
      return new PinIterator(reinterpret_cast<lid_t*>(u(edge)),
                   reinterpret_cast<lid_t*>( toPtr(edge_list[type][id]+1)) );
    }
    return new PinIterator((pin_list[type]+pin_degree_list[type][id]),
                           pin_list[type]+pin_degree_list[type][id+1]);
  }

  void Ngraph::setEdgeWeights(std::vector<wgt_t>& wgts, etype t) {
    assert(!edge_weights[t]);
    if (wgts.size()==0) {
      edge_weights[t] = new wgt_t[num_local_edges[t]];
      for (gid_t i=0;i<num_local_edges[t];i++)
        edge_weights[t][i]=1;
      return;
    }
    assert((lid_t)wgts.size()==num_local_edges[t]);
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

}
