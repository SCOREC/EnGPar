#include "engpar_balancer.h"
#include "engpar_queue.h"
#include <PCU.h>

namespace engpar {

  Balancer::Balancer(agi::Ngraph* g, double f, int v, const char* n) :
    agi::Balancer(g,v,n), factor(f) {
    maxStep =1;
  }
  bool Balancer::runStep(double tolerance) {
    printf("%d Starting runstep\n",PCU_Comm_Self());
    sides = makeSides();
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), sides->print("Sides").c_str());
    vtxWeights = makeVtxWeights(sides);
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), vtxWeights->print("Weights").c_str());
    edgeWeights = new Weights*[graph->numEdgeTypes()];
    for (agi::etype i=0;i<graph->numEdgeTypes();i++) {
      edgeWeights[i] = makeEdgeWeights(sides,i);
    }
    targets = makeTargets(sides,vtxWeights,edgeWeights);
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), targets->print("Targets").c_str());
    Queue* pq = createIterationQueue(graph);
    selector = makeSelector(pq);
    agi::Migration* plan = new agi::Migration;
    wgt_t planW = 0.0;
    for (unsigned int cavSize=2;cavSize<=2;cavSize+=2) {
      planW += selector->select(targets,plan,planW,cavSize);
    }
    graph->migrate(plan);
    return false;
  }
  void Balancer::balance(double tol) {
    if (1 == PCU_Comm_Peers()) return;
    int step = 0;
    while (runStep(tol) && step++<maxStep);
  }
}
