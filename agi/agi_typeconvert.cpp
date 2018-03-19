#include "agi_typeconvert.h"

namespace agi {
  void* toPtr(lid_t t) {
    return (lid_t*)(uintptr_t)(t);
  }

  lid_t fromPtr(GraphVertex* v) {
    lid_t* lid = reinterpret_cast<lid_t*>((char*)(&v));
    return lid[0]-1;
  }
}
