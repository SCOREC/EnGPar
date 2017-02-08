#include "engpar_vtxWeights.h"
#include <PCU.h>
namespace engpar {

  wgt_t getWeight(agi::Ngraph* g) {
    agi::GraphVertex* vtx;
    agi::VertexIterator* vitr = g->begin();
    wgt_t sum = 0.0;
    while ((vtx = g->iterate(vitr))) {
      sum+=g->weight(vtx);
    }
    return sum;
  }

  VtxWeights::VtxWeights(agi::Ngraph* g, Bounds* b) {
    weight = getWeight(g);
    init(b);
  }
  wgt_t VtxWeights::myWeight() {
    return weight;
  }

  void VtxWeights::init(Bounds* b) {
    PCU_Comm_Begin();
    Bounds::iterator itr;
    for (itr=b->begin();itr!=b->end();itr++)
      PCU_COMM_PACK(itr->first,weight);
    PCU_Comm_Send();
    while (PCU_Comm_Listen()) {
      double otherWeight;
      PCU_COMM_UNPACK(otherWeight);
      d[PCU_Comm_Sender()] = otherWeight;
    }
  }

  VtxWeights* makeVtxWeights(agi::Ngraph* g, Bounds* b) {
    return new VtxWeights(g,b);
  }
}
