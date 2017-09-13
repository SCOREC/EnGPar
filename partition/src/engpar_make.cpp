#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"

namespace engpar {

  Sides* makeSides(Input* in) { return new Sides(in);}
  Weights* makeVtxWeights(Input* in, Sides* s) {
    return new Weights(in,s,-1);
  }
  Weights* makeWeights(Input* in, Sides* s, int target) {
    return new Weights(in,s,target);
  }
  Targets* makeTargets(Input* in, Sides* s, Weights* tW,int sT,
                       Weights** cWs, std::vector<wgt_t>& cTs) {
    return new Targets(in,s,tW,sT,cWs,cTs);
  }
}
