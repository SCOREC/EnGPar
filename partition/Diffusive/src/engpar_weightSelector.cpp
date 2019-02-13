#include "engpar_weightSelector.h"
#include <agiMigration.h>
namespace engpar {

  /*Given a graph edge find:
   *  All the non owned vertices surrounding the edge -> others
   *  The vertex that this part owns of the cavity -> <returned value>
   */
  agi::GraphVertex* getWeightCavity(agi::Ngraph* g, agi::GraphEdge* edge, Cavity& others) {
    agi::PinIterator* pitr = g->pins(edge);
    agi::GraphVertex* vtx;
    agi::GraphVertex* v = NULL;
    while ((vtx = g->iterate(pitr))) {
      if (g->owner(vtx)==PCU_Comm_Self()) {
        v = vtx;
      }
      else
        others.insert(vtx);
    }
    g->destroy(pitr);
    return v;
  }

  wgt_t WeightSelector::select(Targets* targets,agi::WeightMigration* plan, wgt_t planW) {

    wgt_t alpha = .1;
    
    q->startIteration();
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      if (planW > targets->total()) break;
      //Create Cavity and peers
      Cavity cav;
      agi::GraphVertex* mine = getWeightCavity(in->g,q->get(itr),cav);
      wgt_t w = in->g->weight(mine) - vSending[mine];
      Cavity::iterator cav_itr;
      for (cav_itr = cav.begin(); cav_itr != cav.end(); ++cav_itr) {
        part_t owner = in->g->owner(*cav_itr);
        if (targets->has(owner) && sending[owner] < targets->get(owner)) {
          int x = w * alpha / cav.size();
          if (x > targets->get(owner) - sending[owner])
            x = targets->get(owner) - sending[owner];
          //Add to plan
          plan->insert(mine, *cav_itr, x);
          sending[owner] += x;
          vSending[mine] += x;
          planW+=x;
        }
      }
      q->addElement(q->get(itr));
    }
    return planW;
  }

  WeightSelector* makeWeightSelector(WeightInput* in, Queue* q) {
    return new WeightSelector(in,q);
  }
}
