#ifndef ENGPAR_VTXWEIGHTS_H
#define ENGPAR_VTXWEIGHTS_H

#include "../engpar.h"
#include "engpar_container.h"
#include "engpar_bounds.h"

namespace engpar {

  class VtxWeights : public Container<wgt_t> {
  public:
    VtxWeights(agi::Ngraph*,Bounds* b);
    wgt_t myWeight();
  private:
    void init(Bounds* b);
    wgt_t weight;
  };

  VtxWeights* makeVtxWeights(agi::Ngraph*, Bounds*);
}

#endif
