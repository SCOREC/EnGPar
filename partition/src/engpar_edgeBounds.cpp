#ifndef ENGPAR_EDGEBOUNDS_H
#define ENGPAR_EDGEBOUNDS_H

#include "engpar_bounds.h"
#include <PCU.h>
namespace {
  class EdgeBounds : public engpar::Bounds {
  public:
    EdgeBounds(agi::Ngraph* g) : engpar::Bounds(g){
      total_boundaries =0;
      init(g);
    }
  private:
    void init(agi::Ngraph* g) {
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g->begin(0);
      while ((edge = g->iterate(eitr))) {
	agi::GraphVertex* pin;
	agi::PinIterator* pitr = g->pins(edge);
	int deg = g->degree(edge);
	for (int i=0;i<deg;i++) {
	  pin = g->iterate(pitr);
	  if (PCU_Comm_Self()!=g->owner(pin)){
	    d[g->owner(pin)]++;
	    total_boundaries++;
	  }
	}
      }
    }
  };
}
namespace engpar {
  Bounds* makeEdgeBounds(agi::Ngraph* graph) {
    return new EdgeBounds(graph);
  }
}

#endif
