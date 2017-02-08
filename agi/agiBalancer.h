#ifndef AGI_BALANCER_H
#define AGI_BALANCER_H

#include "ngraph.h"
namespace agi {

  class Balancer {
  public:
    Balancer(agi::Ngraph* graph_,int verbosity_,const char* name_);

    virtual void balance(double tolerance) =0;
  protected:
    agi::Ngraph* graph;
    int verbosity;
    const char* name;
  };
}
#endif
