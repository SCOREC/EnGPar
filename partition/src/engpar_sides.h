#ifndef ENGPAR_SIDES_H
#define ENGPAR_SIDES_H

#include <ngraph.h>
#include <PCU.h>
#include "engpar_container.h"
#include "../engpar.h"
namespace engpar {
  class Sides : public Container<int>  {
  public:
    Sides(Input* in) {
      agi::Ngraph* graph = in->g;
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = graph->begin(in->sides_edge_type);
      while ((edge = graph->iterate(eitr))) {
        agi::Peers res;
        in->g->getResidence(edge,res);
        if (res.size()>1) {
          agi::Peers::iterator itr;
          for (itr=res.begin();itr!=res.end();itr++) {
            if (*itr!=PCU_Comm_Self()) {
              increment(*itr);
            }
          }
        }
      }
      in->g->destroy(eitr);
    }
  };

  Sides* makeSides(Input* in);
}

#endif
