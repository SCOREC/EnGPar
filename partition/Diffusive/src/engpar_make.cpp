#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"

namespace engpar {

  Sides* makeSides(agi::Ngraph* g, agi::etype t) {return new Sides(g,t);}
  Sides* makeSides(DiffusiveInput* in) { return new Sides(in->g,in->sides_edge_type);}
  Weights* makeVtxWeights(DiffusiveInput* in, Sides* s) {
    return new Weights(in->g,in->countGhosts,s,-1);
  }
  Weights* makeWeights(DiffusiveInput* in, Sides* s, int target) {
    return new Weights(in->g,in->countGhosts,s,target);
  }

  Weights* makeWeights(agi::Ngraph* g, bool countGhosts, Sides* s, int target) {
    return new Weights(g,countGhosts,s,target);
  }
  Targets* makeTargets(double sf, Sides* s, Weights* tW,int sT,
                       Weights** cWs, std::vector<wgt_t>& cTs) {
    return new Targets(sf,s,tW,sT,cWs,cTs);
  }

  Targets* makeTargets(DiffusiveInput* in, Sides* s, Weights* tW,int sT,
                       Weights** cWs, std::vector<wgt_t>& cTs) {
    return new Targets(in->step_factor,s,tW,sT,cWs,cTs);
  }
}
