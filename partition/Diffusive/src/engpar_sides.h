#ifndef ENGPAR_SIDES_H
#define ENGPAR_SIDES_H

#include <ngraph.h>
#include <PCU.h>
#include "engpar_container.h"
#include <engpar_metrics.h>
#include "engpar_diffusive_input.h"

namespace engpar {
  class Sides : public Container<int>  {
  public:
    Sides(agi::Ngraph* g, agi::etype t) {
      agi::Ngraph* graph = g;
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = graph->begin(t);
      while ((edge = graph->iterate(eitr))) {
        agi::Peers res;
        g->getResidence(edge,res);
        if (res.size()>1) {
          agi::Peers::iterator itr;
          for (itr=res.begin();itr!=res.end();itr++) {
            if (*itr!=PCU_Comm_Self()) {
              increment2(*itr);
            }
          }
          my_total++;
        }
      }
      g->destroy(eitr);
    }
  };

  Sides* makeSides(DiffusiveInput* in);
  Sides* makeSides(agi::Ngraph* g, agi::etype t);
  
}

#endif
