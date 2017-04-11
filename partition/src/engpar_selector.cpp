#include "engpar_selector.h"
#include <PCU.h>
namespace engpar {
  
  void getCavity(agi::Ngraph* g, agi::GraphEdge* edge, agi::Migration* plan,
		 Cavity& cav, Peers& peers) {
    agi::PinIterator* pitr = g->pins(edge);
    agi::lid_t deg = g->degree(edge);
    agi::GraphVertex* vtx;
    for (agi::lid_t i =0;i<deg;i++) {
      vtx = g->iterate(pitr);
      if (g->owner(vtx)==PCU_Comm_Self()) {
	if(plan->find(vtx)==plan->end())
	  cav.push_back(vtx);
      }
      else
	peers.push_back(g->owner(vtx));
    }
  }

  wgt_t addCavity(agi::Ngraph* g, Cavity& cav,
		  part_t peer, agi::Migration* plan) {
    Cavity::iterator itr;
    wgt_t w=0.0;
    for (itr = cav.begin();itr!=cav.end();itr++) {
      plan->insert(std::make_pair(*itr,peer));
      w+= g->weight(*itr);
    }
    return w;
  }

  wgt_t Selector::select(Targets* targets, agi::Migration* plan,
			 wgt_t planW, unsigned int cavSize) {
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      //Create Cavity and peers
      Cavity cav;
      Peers peers;
      getCavity(g,*itr,plan,cav,peers);
      //For each peer of cavity
      Peers::iterator itr;
      for (itr = peers.begin();itr!=peers.end();itr++) {
	part_t peer = *itr;
	if (targets->has(peer) &&
	    sending[*itr]<targets->get(peer) &&
	    cav.size()< cavSize) {
	  //  addCavity to plan
	  
	  wgt_t w = addCavity(g,cav,peer,plan);
	  planW+=w;
	  sending[peer]+=w;
	}
      }
    }
    return planW;
  
  }

  void Selector::trim(Targets* targets, agi::Migration* plan) {

  }

  void Selector::cancel(agi::Migration* plan) {

  }
}
