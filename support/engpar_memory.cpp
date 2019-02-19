#include "engpar_support.h"

#ifdef __bgq__
#include <spi/include/kernel/memory.h>

double EnGpar_Peak_Memory()
{
  uint64_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  return heap;
}

#elif defined (__linux__)
#include <malloc.h>
double EnGPar_Peak_Memory()
{
  return mallinfo().arena;
}

#else

double EnGPar_Peak_Memory()
{
  printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}

#endif
