#ifndef ENGPAR_H
#define ENGPAR_H
#include <agiBalancer.h>
#include <ngraph.h>
#include <vector>
#include "engpar_input.h"

namespace engpar {
  typedef agi::wgt_t wgt_t;
  typedef agi::part_t part_t;
  typedef std::vector<agi::GraphEdge*> Queue;
  
  wgt_t getWeight(agi::Ngraph*,int);
  //Balancers:
  agi::Balancer* makeVtxBalancer(agi::Ngraph*, double stepFactor=0.1,
                                        int verbosity=0);
  agi::Balancer* makeBalancer(Input*,int verbosity=0);

  double EnGPar_Get_Imbalance(wgt_t);
  //Partition info:
  void evaluatePartition(agi::Ngraph*);
}
#endif
