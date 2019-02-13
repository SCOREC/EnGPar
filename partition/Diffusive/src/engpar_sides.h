#ifndef ENGPAR_SIDES_H
#define ENGPAR_SIDES_H

#include <ngraph.h>
#include <PCU.h>
#include "engpar_container.h"
#include <engpar_metrics.h>
#include "engpar_diffusive_input.h"

namespace engpar {
  class Sides : public Container<int>  {
  public:
    Sides(agi::Ngraph* g, agi::etype t) {
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g->begin(t);
      //For each edge
      while ((edge = g->iterate(eitr))) {
        agi::Peers res;
        g->getResidence(edge,res);
	//Look at edges that are cut
        if (res.size() > 1) {
          agi::Peers::iterator itr;
          for (itr=res.begin();itr!=res.end();itr++) {
            if (*itr != PCU_Comm_Self()) {
	      //Add the weight of the edge for the side of the neighbor
              (*this)[*itr]+=g->weight(edge);
            }
          }
	  //Total the weight of all cut edges
          my_total+= g->weight(edge);
        }
      }
      g->destroy(eitr);
    }
  };

  Sides* makeSides(DiffusiveInput* in);
  Sides* makeSides(agi::Ngraph* g, agi::etype t);
  
}

#endif
