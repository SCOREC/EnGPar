#include "agi_typeconvert.h"

namespace agi {
  void* toPtr(lid_t t) {
    return (lid_t*)(uintptr_t)(t);
  }

  lid_t fromPtr(GraphVertex* v) {
    lid_t* lid = reinterpret_cast<lid_t*>((char*)(&v));
#ifndef ENGPAR_BIG_ENDIAN
    return lid[0]-1;
#elif LOCAL_ID_TYPE == int64_t
    return lid[0]-1;
#else
    return lid[1]-1;
#endif
  }
}
