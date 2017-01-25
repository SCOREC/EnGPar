#ifndef __EDGE_ITERATOR__H__
#define __EDGE_ITERATOR__H__
#include "agi.h"

namespace agi {

class Ngraph;

class EdgeIterator {
  friend class Ngraph;
 private:
  gid_t* loc;
  gid_t* end;
 EdgeIterator(etype t, int nt, lid_t* l,lid_t deg) : loc((gid_t*)(1+t+nt*(uintptr_t)l)),
    end((gid_t*)(1+t+nt*((uintptr_t)l+deg))){}
};
  
}


#endif
