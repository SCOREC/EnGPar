#include "agi_typeconvert.h"

namespace agi {
  void* toPtr(lid_t t) {
    return (lid_t*)(uintptr_t)(t);
  }

  lid_t fromPtr(GraphVertex* v) {
    return reinterpret_cast<lid_t*>((char*)(&v))[0]-1;
  }
}
