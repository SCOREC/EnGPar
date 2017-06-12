#include "engpar_input.h"

namespace engpar {
  
Input::Input(agi::Ngraph* g_) {
  g=g_;
  maxIterations = 1000;
  maxIterationsPerType = 100;
  step_factor=0.1;
  sides_edge_type = 0;
  selection_edge_type = 0;
  countGhosts =false;
}

}
