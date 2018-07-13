#ifndef __COLORING_INPUT_H__
#define __COLORING_INPUT_H__

#include "../engpar_input.h"
#include <ngraph.h>
#include <PCU.h>

namespace engpar {
  class ColoringInput : public Input {
  public:
  //  ColoringInput(agi::Ngraph* g_, agi::lid_t t) : Input(g_) {
  //    assert( (t == VTX_TYPE || t < g->numEdgeTypes()) &&
  //            "invalid primary type for coloring input");
  //    primaryType = t;
  //    if (t==VTX_TYPE) {
  //      edgeType = 0;
  //    }
  //    else
  //      edgeType = t; 
  //  }
    ColoringInput(agi::Ngraph* g_, agi::lid_t t, bool vtx) : Input(g_) {
      assert( (t < g->numEdgeTypes()) &&
              "invalid primary type for coloring input");
      if (vtx)
        primaryType = VTX_TYPE;
      else
        primaryType = t;
      edgeType = t;
    }
    agi::lid_t primaryType;
    agi::lid_t edgeType;
  };
}

#endif
