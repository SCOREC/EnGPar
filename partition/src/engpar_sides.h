#ifndef ENGPAR_SIDES_H
#define ENGPAR_SIDES_H

#include <ngraph.h>
#include "engpar_container.h"
namespace engpar {
  class Sides : public Container<int>  {
  public:
    Sides() {}
  };

  //  Sides* makeEdgeSides(agi::Ngraph*);
}

#endif
