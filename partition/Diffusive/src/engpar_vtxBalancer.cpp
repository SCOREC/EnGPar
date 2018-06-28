#include <engpar_support.h>
#include <engpar_metrics.h>
#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_balancer.h"
#include "engpar_targets.h"
#include <PCU.h>
namespace engpar {
  class VtxBalancer : public engpar::Balancer {
  public:
    VtxBalancer(agi::Ngraph*& g, double f, int v, double tol)
      : Balancer(g,f,v,"Vtx") {
      DiffusiveInput* inp = dynamic_cast<DiffusiveInput*>(input);
      inp->addPriority(-1,tol);
    }
    ~VtxBalancer() {}
  };


  void balanceVertices(agi::Ngraph*& g, double tol, double f, int v) {
    Balancer* balancer = new VtxBalancer(g,f,v,tol);
    balancer->balance();
    delete balancer;
  }
}
