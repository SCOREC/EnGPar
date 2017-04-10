#ifndef ENGPAR_BALANCER_H
#define ENGPAR_BALANCER_H

#include "../engpar.h"
#include <agiBalancer.h>
#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"
#include "engpar_selector.h"

namespace engpar {

  class Balancer : public agi::Balancer{
  public:
    Balancer(agi::Ngraph* graph_, double factor_, int verbosity_,
	     const char* name_);
    virtual ~Balancer() {}
    virtual bool runStep(double tolerance);
    void balance(double tolerance);
    virtual Sides* makeSides()=0;
    virtual Weights* makeVtxWeights(Sides* s)=0;
    virtual Weights* makeEdgeWeights(Sides* s, agi::etype i)=0;
    virtual Targets* makeTargets(Sides* s, Weights* vtxW,
				 Weights** edgeW)=0;
    virtual Selector* makeSelector(Queue*)=0;
  protected:
    double factor;
    int maxStep;
    Sides* sides;
    Weights* vtxWeights,**edgeWeights;
    Targets* targets;
    Selector* selector;
    
  };
}

#endif
