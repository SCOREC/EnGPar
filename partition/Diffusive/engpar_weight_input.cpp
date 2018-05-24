#include "engpar_weight_input.h"
namespace engpar {

  WeightInput::WeightInput(agi::Ngraph* g_,double t, double sf, agi::etype et) : Input(g_) {
    maxIterations = 100;
    tol = t;
    step_factor=sf;
    primary_edge_type = et;
  }
}
