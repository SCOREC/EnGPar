#include "ngraph.h"
#include <PCU.h>
namespace agi {

  bool checkValidity(Ngraph* g) {
    lid_t i =0;
    VertexIterator* vitr = g->begin();
    GraphVertex* vtx;
    while ((vtx = g->iterate(vitr))) {
      i++;
    }
    if(g->numLocalVtxs()!=i) return false;
    lid_t ig = PCU_Add_Long(i);
    if(g->numGlobalVtxs()!=ig) return false;
    vitr = g->begin();
    while ((vtx = g->iterate(vitr))) {
      for (etype t = 0;t<g->numEdgeTypes();t++) {
        GraphEdge* edge;
        EdgeIterator* eitr = g->edges(vtx,t);
        lid_t j=0;
        while ((edge = g->iterate(eitr))) {
          j++;
        }
        g->destroy(eitr);
        if(g->degree(vtx,t)!=j) return false;
      }
    }

    vitr = g->begin();
    while ((vtx = g->iterate(vitr))) {
      for (etype t = 0;t<g->numEdgeTypes();t++) {
        GraphEdge* edge;
        EdgeIterator* eitr = g->edges(vtx,t);
        while ((edge = g->iterate(eitr))) {
          PinIterator* pitr = g->pins(edge);
          GraphVertex* other;
          while ((other = g->iterate(pitr))) {
            if (g->isEqual(other,vtx))
              break;
          }
          if(!other) return false;
          g->destroy(pitr);
        }
        g->destroy(eitr);
      }
    }

    for (etype t=0;t<g->numEdgeTypes();t++) {
      EdgeIterator* eitr = g->begin(t);
      GraphEdge* edge;
      i=0;
      lid_t j=0;
      lid_t k=0;
      lid_t l=0;
      while ((edge = g->iterate(eitr))) {
        i++;
        l+=g->degree(edge);
        Peers res;
        g->getResidence(edge,res);
        Peers::iterator itr;
        part_t owner=PCU_Comm_Self();
        for (itr=res.begin();itr!=res.end();itr++) {
          if (*itr<owner)
            break;
        }
        if (itr==res.end()||!g->isHyper()){
          j++;
          k+=g->degree(edge);
        }
      }
      g->destroy(eitr);
      if(i!=g->numLocalEdges(t)) return false;
      j = PCU_Add_Long(j);
      if(j!=g->numGlobalEdges(t)) return false;
      k = PCU_Add_Long(k);
      if(k!=g->numGlobalPins(t)) return false;
      if(l!=g->numLocalPins(t)) return false;
      eitr = g->begin(t);
      while ((edge = g->iterate(eitr))) {
        GraphVertex* other;
        PinIterator* pitr = g->pins(edge);
        j=0;
        while ((other = g->iterate(pitr))) {
          j++;
        }
        g->destroy(pitr);
        if(g->degree(edge)!=j) return false;
      }
      g->destroy(eitr);
    }
    //Ensure all ghost vertices have valid owners
    GhostIterator* g_itr = g->beginGhosts();
    GraphVertex* gv;
    while ((gv = g->iterate(g_itr))) {
      if(g->owner(gv)==PCU_Comm_Self()) return false;
      if(g->owner(gv)>=PCU_Comm_Peers()) return false;
    }
    
    return true;
  }
  
  
}
