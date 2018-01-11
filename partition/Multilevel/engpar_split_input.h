#ifndef __SPLIT_INPUT_H__
#define __SPLIT_INPUT_H__

#include "../engpar_input.h"
#include <ngraph.h>

namespace engpar {
  class SplitInput : public Input {
  public:
    SplitInput(agi::Ngraph*& g_) : Input(g_) {
      smallComm = largeComm = MPI_COMM_WORLD;
      isOriginal=true;
      split_factor = 1;
      tolerance = 1.05;
      edge_type = 0;
    }
    MPI_Comm smallComm;
    MPI_Comm largeComm;
    bool isOriginal;
    int split_factor;
    double tolerance;
    agi::etype edge_type;
    
  };
}

#endif
