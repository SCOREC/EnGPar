#ifndef ENGPAR_TARGETS_H
#define ENGPAR_TARGETS_H

#include "../engpar.h"
#include "engpar_container.h"
#include "engpar_bounds.h"
#include "engpar_vtxWeights.h"
namespace engpar {

  class Targets : public Container<wgt_t> {
  public:
    virtual wgt_t total()=0;
  };

  Targets* makeTargets(Bounds* b, VtxWeights* w,double factor);
}

#endif
