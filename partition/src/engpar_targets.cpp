#include "engpar_targets.h"

namespace {
  class WeightTargets : public engpar::Targets {
  public:
    WeightTargets(engpar::Bounds* b, engpar::VtxWeights* w, double factor) {
      init(b,w,factor);
    }
    engpar::wgt_t total() {
      return totW;
    }
  private:
    engpar::wgt_t totW;
    void init(engpar::Bounds* b, engpar::VtxWeights* w, double factor) {
      totW=0.0;
      engpar::Bounds::iterator itr;
      for (itr = b->begin();itr!=b->end();itr++) {
	int neighbor = itr->first;
	engpar::wgt_t myW = w->myWeight();
	engpar::wgt_t neighborW = w->get(neighbor);
	if (myW>neighborW) {
	  engpar::wgt_t diff = myW-neighborW;
	  engpar::wgt_t sideFraction = itr->second;
	  sideFraction /= b->total();
	  engpar::wgt_t scaledW = diff * sideFraction* factor;
	  set(neighbor,scaledW);
	  totW+=scaledW;
	}
      }
    }
  };
}
namespace engpar {
  Targets* makeTargets(Bounds* b, VtxWeights* w, double factor) {
    return new WeightTargets(b,w,factor);
  }
}
