#ifndef ENGPAR_SIDES_H
#define ENGPAR_SIDES_H

#include <ngraph.h>
#include <PCU.h>
#include "engpar_container.h"
#include "engpar_input.h"
namespace engpar {
  class Sides : public Container<int>  {
  public:
    Sides(Input* in) {
      agi::Ngraph* graph = in->g;
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = graph->begin(in->sides_edge_type);
      while ((edge = graph->iterate(eitr))) {
        agi::GraphVertex* pin;
        agi::PinIterator* pitr = graph->pins(edge);
        int deg = graph->degree(edge);
        for (int i=0;i<deg;i++) {
          pin = graph->iterate(pitr);
          agi::part_t owner = graph->owner(pin);
          if (PCU_Comm_Self()!=owner){
            increment(owner);
          }
        }
      }
    }
  };

  Sides* makeSides(Input* in);
}

#endif
