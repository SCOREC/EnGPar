#ifndef ENGPAR_H
#define ENGPAR_H
#include <agiBalancer.h>
#include <ngraph.h>
#include <vector>
namespace engpar {
  typedef agi::wgt_t wgt_t;
  typedef agi::part_t part_t;
  typedef std::vector<agi::GraphEdge*> Queue;

  //Balancers:
  agi::Balancer* makeVtxBalancer(agi::Ngraph* g, double stepFactor=0.1,
					int verbosity=0);


  //Partition info:
  void evaluatePartition(agi::Ngraph* g);
}
#endif
