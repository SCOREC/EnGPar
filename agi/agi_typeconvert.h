#ifndef AGI_TYPE_CONVERT_H_
#define AGI_TYPE_CONVERT_H_

#include "agi.h"

/**
 * use these with extreme care...
 **/

namespace agi {
  void* toPtr(lid_t t);
  lid_t fromPtr(GraphVertex* v);
}

#endif
