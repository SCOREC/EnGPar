#include "../engpar.h"
#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_balancer.h"
#include "engpar_targets.h"
#include <unordered_set>
#include <PCU.h>
namespace engpar {
  class VtxBalancer : public engpar::Balancer {
  public:
    VtxBalancer(agi::Ngraph* g, double f, int v)
      : Balancer(g,f,v,"Vtx") {
      input->priorities.push_back(-1);
    }
    ~VtxBalancer() {}
  };
}

namespace engpar {
  agi::Balancer* makeVtxBalancer(agi::Ngraph* g, double f, int v) {
    return new VtxBalancer(g,f,v);
  }
}
