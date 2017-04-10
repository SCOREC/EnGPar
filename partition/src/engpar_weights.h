#ifndef ENGPAR_WEIGHTS_H
#define ENGPAR_WEIGHTS_H

#include "../engpar.h"
#include "engpar_container.h"
#include "engpar_sides.h"

namespace engpar {
  class Balancer;
  
  class Weights : public Container<wgt_t> {
  public:
    const wgt_t& myWeight() const {return my_weight;}
    void addWeight(wgt_t w) {my_weight+=w;}
  private:
    wgt_t my_weight;
  };

  //VtxWeights* makeVtxWeights(agi::Ngraph*, Bounds*);
}

#endif
