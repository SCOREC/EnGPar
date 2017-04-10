#include "engpar_queue.h"
#include <PCU.h>

namespace engpar {
  Queue* createIterationQueue(agi::Ngraph* g) {
    Queue* q = new Queue;
    agi::GraphEdge* edge;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((edge = g->iterate(eitr))) {
      agi::PinIterator* pitr = g->pins(edge);
      agi::lid_t degree = g->degree(edge);
      for (agi::lid_t i=0;i<degree;i++) {
	agi::GraphVertex* pin = g->iterate(pitr);
	if (g->owner(pin)!=PCU_Comm_Self()) {
	  q->push_back(edge);
	  break;
	}
      }
    }
    return q;
  }
}
