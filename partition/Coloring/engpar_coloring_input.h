#ifndef __COLORING_INPUT_H__
#define __COLORING_INPUT_H__

#include "../engpar_input.h"
#include <ngraph.h>
#include <PCU.h>

namespace engpar {
  class ColoringInput : public Input {
  public:
    ColoringInput(agi::Ngraph* g_) : Input(g_) {
      primaryType = 0;
    }
    ColoringInput(agi::Ngraph* g_, agi::lid_t t) : Input(g_) {
      assert( (t == -1 || t < g->numEdgeTypes()) &&
              "invalid primary type for coloring input");
      primaryType = t;
    }
    agi::lid_t primaryType;
  };
}

#endif
