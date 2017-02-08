#ifndef ENGPAR_BOUNDS_H
#define ENGPAR_BOUNDS_H

#include <ngraph.h>
#include "engpar_container.h"
namespace engpar {
  class Bounds : public Container<int>  {
  public:
    Bounds(agi::Ngraph* g) {}
    int total() {return total_boundaries;}
  protected:
    int total_boundaries;    
  };

  Bounds* makeEdgeBounds(agi::Ngraph*);
}

#endif
