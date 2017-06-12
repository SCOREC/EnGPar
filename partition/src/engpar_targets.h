#ifndef ENGPAR_TARGETS_H
#define ENGPAR_TARGETS_H

#include "../engpar.h"
#include "engpar_container.h"
#include "engpar_weights.h"
namespace engpar {

  class Targets : public Container<wgt_t> {
  public:
    Targets(Input* in, Sides* s, Weights* vtxW,
	    Weights** edgeWs) {
      Sides::iterator itr;
      for (itr = s->begin();itr!=s->end();itr++) {
	int neighbor = itr->first;
	engpar::wgt_t myW = vtxW->myWeight();
	engpar::wgt_t neighborW = vtxW->get(neighbor);
	if (myW>neighborW) {
	  engpar::wgt_t diff = myW-neighborW;
	  engpar::wgt_t sideFraction = itr->second;
	  sideFraction /= s->total();
	  engpar::wgt_t scaledW = diff * sideFraction* in->step_factor;
	  set(neighbor,scaledW);
	}
      }
    }
  };
  Targets* makeTargets(Input* in, Sides* s, Weights* vW,Weights** eWs);
}

#endif
