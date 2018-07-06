#include "agi_typeconvert.h"

#define ENGPAR_INT_IF_int32_t 1
#define ENGPAR_INT_IF_int64_t 0
                               
#define ENGPAR_COMBINE(a,b) a ## b
#define ENGPAR_INT_IF(val) ENGPAR_COMBINE(ENGPAR_INT_IF_, val)

namespace agi {
  void* toPtr(lid_t t) {
    return (lid_t*)(uintptr_t)(t);
  }

  lid_t fromPtr(GraphVertex* v) {
    lid_t* lid = reinterpret_cast<lid_t*>((char*)(&v));
#ifndef ENGPAR_BIG_ENDIAN
    return lid[0]-1;
#elif ENGPAR_INT_IF(LOCAL_ID_TYPE) == 0
    return lid[0]-1;
#else
    return lid[1]-1;
#endif
  }
}
