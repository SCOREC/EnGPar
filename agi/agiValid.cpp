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
    assert(g->numLocalVtxs()==i);
    lid_t ig = PCU_Add_Long(i);
    assert(g->numGlobalVtxs()==ig);
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
        assert(g->degree(vtx,t)==j);
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
          assert(other);
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
      assert(i==g->numLocalEdges(t));
      j = PCU_Add_Long(j);
      assert(j==g->numGlobalEdges(t));
      k = PCU_Add_Long(k);
      assert(k==g->numGlobalPins(t));
      assert(l==g->numLocalPins(t));
      eitr = g->begin(t);
      while ((edge = g->iterate(eitr))) {
        GraphVertex* other;
        PinIterator* pitr = g->pins(edge);
        j=0;
        while ((other = g->iterate(pitr))) {
          j++;
        }
        g->destroy(pitr);
        assert(g->degree(edge)==j);
      }
      g->destroy(eitr);
    }
    //Ensure all ghost vertices have valid owners
    GhostIterator* g_itr = g->beginGhosts();
    GraphVertex* gv;
    while ((gv = g->iterate(g_itr))) {
      assert(g->owner(gv)!=PCU_Comm_Self());
      assert(g->owner(gv)<PCU_Comm_Peers());
    }
    
    return true;
  }
  
  
}
