#ifndef ENGPAR_H
#define ENGPAR_H
#include <agiBalancer.h>
#include <ngraph.h>

namespace engpar {
typedef agi::wgt_t wgt_t;
typedef agi::part_t part_t;


  agi::Balancer* makeVtxBalancer(agi::Ngraph* g, double stepFactor=0.1,
					int verbosity=0);
}
#endif
