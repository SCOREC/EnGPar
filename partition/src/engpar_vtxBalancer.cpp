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
      
    }
    ~VtxBalancer() {}

    //Count the number of pins across part boundaries
    virtual Sides* makeSides() {
      Sides* s = new Sides;
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = graph->begin(0);
      while ((edge = graph->iterate(eitr))) {
	agi::GraphVertex* pin;
	agi::PinIterator* pitr = graph->pins(edge);
	int deg = graph->degree(edge);
	for (int i=0;i<deg;i++) {
	  pin = graph->iterate(pitr);
	  part_t owner = graph->owner(pin);
	  if (PCU_Comm_Self()!=owner){
	    s->increment(owner);
	  }
	}
      }
      return s;
    }
    
    virtual Weights* makeVtxWeights(Sides* s) {
      Weights* w = new Weights();
      //calculate the total weight of the vertices
      w->addWeight(getWeight(graph,-1));

      //Share weight with all neighbors
      PCU_Comm_Begin();
      Sides::iterator itr;
      for (itr=s->begin();itr!=s->end();itr++) 
	PCU_COMM_PACK(itr->first,w->myWeight());
      
      PCU_Comm_Send();
      while (PCU_Comm_Listen()) {
	double otherWeight;
	PCU_COMM_UNPACK(otherWeight);
	w->set(PCU_Comm_Sender(),otherWeight);
      }
      return w;

    }
    //Only balancing vertices so we won't have edge weights
    virtual Weights* makeEdgeWeights(Sides* s, agi::etype i) {
      return NULL;
    }

    //Calculates the amount fo weight to send to each neighbor
    virtual Targets* makeTargets(Sides* s, Weights* vtxW,
				 Weights** edgeW) {
      Targets* t = new Targets();
      Sides::iterator itr;
      for (itr = s->begin();itr!=s->end();itr++) {
	int neighbor = itr->first;
	engpar::wgt_t myW = vtxW->myWeight();
	engpar::wgt_t neighborW = vtxW->get(neighbor);
	if (myW>neighborW) {
	  engpar::wgt_t diff = myW-neighborW;
	  engpar::wgt_t sideFraction = itr->second;
	  sideFraction /= s->total();
	  engpar::wgt_t scaledW = diff * sideFraction* factor;
	  t->set(neighbor,scaledW);
	}
      }
      return t;
    }
    virtual Selector* makeSelector(Queue* q) {
      return new Selector(graph,q,&completed_dimensions,&completed_weights);
    }

  };
}

namespace engpar {
  agi::Balancer* makeVtxBalancer(agi::Ngraph* g, double f, int v) {
    return new VtxBalancer(g,f,v);
  }
}
