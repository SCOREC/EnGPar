#include "ngraph.h"
#include <PCU.h>
#include <engpar_support.h>
namespace agi {

  bool checkValidity(Ngraph* g) {
    lid_t i =0;
    VertexIterator* vitr = g->begin();
    GraphVertex* vtx;
    while ((vtx = g->iterate(vitr))) {
      i++;
    }
    if(g->numLocalVtxs()!=i) {
      EnGPar_Error_Message("Number of local vertices does not match the number of "
                           "vertices on process %d\n", PCU_Comm_Self());
      return false;
    }
    lid_t ig = PCU_Add_Long(i);
    if(g->numGlobalVtxs()!=ig) {
      EnGPar_Error_Message("Number of global vertices does not match the sum of "
                           "the number of vertices on each process\n");
      return false;
    }
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
        if(g->degree(vtx,t)!=j) {
          EnGPar_Error_Message("Process: %d Degree of vertex %ld does not match the "
                               "number of edges adjacent to it.", PCU_Comm_Self(),vtx);
          return false;
        }
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
          if(!other) {
            EnGPar_Error_Message("Process %d: Vertex, %ld, is connected to edge, %ld, "
                                 "but the edge is not connected to the vertex\n",
                                 PCU_Comm_Self(),g->globalID(vtx),g->globalID(edge));
            return false;
          }
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
      if(i!=g->numLocalEdges(t)) {
        EnGPar_Error_Message("Process %d: Number of local edges does not match "
                             "the number of actual edges\n",PCU_Comm_Self());
        return false;
      }
      j = PCU_Add_Long(j);
      if(j!=g->numGlobalEdges(t)) {
        EnGPar_Error_Message("Number of global edges does not match the sum of "
                             "the number of edges on each process\n");
        return false;
      }
      if(l!=g->numLocalPins(t)) {
        EnGPar_Error_Message("Process %d: Number of local pins does not match "
                             "the number of actual pins\n",PCU_Comm_Self());
        return false;
      }
      k = PCU_Add_Long(k);
      if(k!=g->numGlobalPins(t)) {
        EnGPar_Error_Message("Number of global pins does not match the sum of "
                             "the number of pins on each process\n");
        return false;
      }
      eitr = g->begin(t);
      while ((edge = g->iterate(eitr))) {
        GraphVertex* other;
        PinIterator* pitr = g->pins(edge);
        j=0;
        while ((other = g->iterate(pitr))) {
          j++;
        }
        g->destroy(pitr);
        if(g->degree(edge)!=j) {
          EnGPar_Error_Message("Process: %d Degree of edge %ld does not match the "
                               "number of vertices that it connects.", PCU_Comm_Self(),vtx);
          return false;
        }
      }
      g->destroy(eitr);
    }
    //Ensure all ghost vertices have valid owners
    GhostIterator* g_itr = g->beginGhosts();
    GraphVertex* gv;
    while ((gv = g->iterate(g_itr))) {
      if(g->owner(gv)==PCU_Comm_Self()) {
        EnGPar_Error_Message("Process %d: invalid ghost vertex, %d, that is owned by itself\n",
                             PCU_Comm_Self(),g->globalID(gv));
        return false;
      }
      if(g->owner(gv)>=PCU_Comm_Peers()) {
        EnGPar_Error_Message("Process %d: ghost vertex, %d, is owned by a process that does "
                             "not exist\n", PCU_Comm_Self(),g->globalID(gv));
                             
        return false;
      }
    }
    
    return true;
  }
  
  
}
