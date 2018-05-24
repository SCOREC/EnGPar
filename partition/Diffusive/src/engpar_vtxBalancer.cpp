#include <engpar_support.h>
#include "../engpar.h"
#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_balancer.h"
#include "engpar_targets.h"
#include <PCU.h>
namespace engpar {
  class VtxBalancer : public engpar::Balancer {
  public:
    VtxBalancer(agi::Ngraph*& g, double f, int v)
      : Balancer(g,f,v,"Vtx") {
      DiffusiveInput* inp = dynamic_cast<DiffusiveInput*>(input);
      inp->addPriority(-1,1.1);
    }
    ~VtxBalancer() {}
  };


  void balanceVertices(agi::Ngraph*& g, double tol, double f, int v) {
    Balancer* balancer = new VtxBalancer(g,f,v);
    balancer->balance(tol);
    delete balancer;
  }
  /*
  agi::Balancer* makeVtxBalancer(agi::Ngraph*& g, double f, int v) {
    if (EnGPar_Is_Log_Open()) {
      char message[25];
      sprintf(message,"makeBalancer\n");
      EnGPar_Log_Function(message);
      EnGPar_End_Function();
    }
    return new VtxBalancer(g,f,v);
  }
  */
}
