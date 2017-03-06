#ifndef ENGPAR_BALANCER_H
#define ENGPAR_BALANCER_H

#include "../engpar.h"
#include <agiBalancer.h>
namespace engpar {

  class Balancer : public agi::Balancer{
  public:
    Balancer(agi::Ngraph* graph_, double factor_, int verbosity_,
	     const char* name_);
    virtual ~Balancer() {}
    virtual bool runStep(double tolerance)=0;
    void balance(double tolerance);
  protected:
    double factor;
    int maxStep;
  };
}

#endif
