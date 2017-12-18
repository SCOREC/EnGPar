#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"

namespace engpar {

  Sides* makeSides(DiffusiveInput* in) { return new Sides(in);}
  Weights* makeVtxWeights(DiffusiveInput* in, Sides* s) {
    return new Weights(in,s,-1);
  }
  Weights* makeWeights(DiffusiveInput* in, Sides* s, int target) {
    return new Weights(in,s,target);
  }
  Targets* makeTargets(DiffusiveInput* in, Sides* s, Weights* tW,int sT,
                       Weights** cWs, std::vector<wgt_t>& cTs) {
    return new Targets(in,s,tW,sT,cWs,cTs);
  }
}
