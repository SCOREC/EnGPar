#ifndef ENGPAR_WEIGHT_SELECTOR_H__
#define ENGPAR_WEIGHT_SELECTOR_H__

#include <engpar_metrics.h>
#include "engpar_sides.h"
#include "engpar_targets.h"
#include "engpar_queue.h"
#include "engpar_weight_input.h"
#include "engpar_selector.h"

namespace agi {
  class WeightMigration;
}
namespace engpar {

  typedef std::unordered_map<part_t,wgt_t> Sending;
  typedef std::unordered_map<agi::GraphVertex*,wgt_t> VertexSending;
  class WeightSelector {

  public:
  WeightSelector(WeightInput* in_, Queue* queue) : in(in_), q(queue) {}

    wgt_t select(Targets* targets,agi::WeightMigration* plan, wgt_t planW);

  private:

    WeightInput* in;
    Queue* q;
    Sending sending;
    VertexSending vSending;
  };

  WeightSelector* makeWeightSelector(WeightInput* in, Queue* q);
}

#endif
