#include "agiBalancer.h"

namespace agi {

  Balancer::Balancer(agi::Ngraph* g,int v,const char* n)
    : graph(g), verbosity(v), name(n) {
  }

}
