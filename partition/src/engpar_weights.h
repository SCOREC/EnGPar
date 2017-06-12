#ifndef ENGPAR_WEIGHTS_H
#define ENGPAR_WEIGHTS_H

#include "../engpar.h"
#include "engpar_container.h"
#include "engpar_sides.h"
#include "engpar_input.h"

namespace engpar {
  class Balancer;
  
  class Weights : public Container<wgt_t> {
  public:
    Weights(Input* in,Sides* s, int target_dimension) {
      //calculate the total weight of the vertices
      my_weight = getWeight(in->g,target_dimension);
      //Share weight with all neighbors
      PCU_Comm_Begin();
      Sides::iterator itr;
      for (itr=s->begin();itr!=s->end();itr++) 
	PCU_COMM_PACK(itr->first,myWeight());
      PCU_Comm_Send();
      while (PCU_Comm_Listen()) {
	double otherWeight;
	PCU_COMM_UNPACK(otherWeight);
	set(PCU_Comm_Sender(),otherWeight);
      }
    }
    const wgt_t& myWeight() const {return my_weight;}
    void addWeight(wgt_t w) {my_weight+=w;}
  private:
    wgt_t my_weight;
  };

  Weights* makeVtxWeights(Input* in, Sides* s);
}

#endif
