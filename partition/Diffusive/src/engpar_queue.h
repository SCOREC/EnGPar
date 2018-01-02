#include <ngraph.h>
#include "../engpar.h"

namespace engpar {
  Queue* createIterationQueue(agi::Ngraph*);
  Queue* createDistanceQueue(DiffusiveInput*);
}
