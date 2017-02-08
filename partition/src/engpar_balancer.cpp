#include "engpar_balancer.h"
#include <PCU.h>
namespace engpar {

  Balancer::Balancer(agi::Ngraph* g, double f, int v, const char* n) :
    agi::Balancer(g,v,n), factor(f) {
    maxStep =1;
  }

  void Balancer::balance(double tol) {
    if (1 == PCU_Comm_Peers()) return;
    int step = 0;
    while (runStep(tol) && step++<maxStep);
  }
}
