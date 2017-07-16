#include <ngraph.h>
#include "../engpar.h"

namespace engpar {
  Queue* createIterationQueue(agi::Ngraph* g);
  Queue* createDistanceQueue(agi::Ngraph* g);
}
