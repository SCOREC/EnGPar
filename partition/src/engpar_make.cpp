#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"

namespace engpar {

  Sides* makeSides(Input* in) { return new Sides(in);}
  Weights* makeVtxWeights(Input* in, Sides* s) {
    return new Weights(in,s,-1);
  }
  Targets* makeTargets(Input* in, Sides* s, Weights* vW,Weights** eWs) {
    return new Targets(in,s,vW,eWs);
  }
}
