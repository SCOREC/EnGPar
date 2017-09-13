#ifndef ENGPAR_TARGETS_H
#define ENGPAR_TARGETS_H

#include "../engpar.h"
#include "engpar_container.h"
#include "engpar_weights.h"
namespace engpar {

  class Targets : public Container<wgt_t> {
  public:
    Targets(Input* in, Sides* s, Weights* targetW, int sideTol,
            Weights** completedWs,std::vector<wgt_t>& completedTolerances) {
      Sides::iterator itr;
      for (itr = s->begin();itr!=s->end();itr++) {
        int neighbor = itr->first;
        bool canSend=s->get(neighbor)<sideTol;
        for (unsigned int i=0;i<completedTolerances.size();i++) {
          if (completedWs[i]->get(neighbor)>=completedTolerances[i])
            canSend=false;
        }
        if (!canSend)
          continue;
        wgt_t myW = targetW->myWeight();
        wgt_t neighborW = targetW->get(neighbor);
        if (myW>neighborW) {
          wgt_t diff = myW-neighborW;
          wgt_t sideFraction = itr->second;
          sideFraction /= s->total();
          wgt_t scaledW = diff * sideFraction* in->step_factor;
          set(neighbor,scaledW);
        }
      }
    }
  };
  Targets* makeTargets(Input* in, Sides* s, Weights* tW,int sT,
                       Weights** cWs,std::vector<wgt_t>& cTs);
}

#endif
