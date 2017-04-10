#ifndef ENGPAR_SELECTOR_H
#define ENGPAR_SELECTOR_H

#include "../engpar.h"
#include "engpar_sides.h"
#include "engpar_targets.h"
#include "engpar_queue.h"

namespace engpar {

  typedef std::vector<agi::GraphVertex*> Cavity;
  typedef std::vector<part_t> Peers;
  typedef std::unordered_map<part_t,wgt_t> Sending;
  
  class Selector {
  public:
    Selector(agi::Ngraph* graph, Queue* queue) : g(graph), q(queue) {}

    virtual wgt_t select(Targets* targets,agi::Migration* plan,
			wgt_t planW, unsigned int cavSize);
    virtual void trim(Targets* targets, agi::Migration* plan);
    virtual void cancel(agi::Migration* plan);
  protected:
    
    agi::Ngraph* g;
    Sending sending;
    Queue* q;
  };
}

#endif
