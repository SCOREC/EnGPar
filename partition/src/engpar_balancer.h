#ifndef ENGPAR_BALANCER_H
#define ENGPAR_BALANCER_H

#include "../engpar.h"
#include <agiBalancer.h>
#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"
#include "engpar_selector.h"

namespace engpar {
  wgt_t getMaxWeight(agi::Ngraph*, int);
  wgt_t getAvgWeight(agi::Ngraph*, int);
  class Balancer : public agi::Balancer{
  public:
    Balancer(agi::Ngraph* graph_, double factor_, int verbosity_,
	     const char* name_);
    Balancer(Input* input_,int verbosity_,const char* name_);
    virtual ~Balancer() {}
    virtual bool runStep(double tolerance);
    void balance(double tolerance);
    
  protected:
    int target_dimension;
    std::vector<int> completed_dimensions;
    std::vector<wgt_t> completed_weights;
    Input* input;
    double times[2];
  };
}

#endif
