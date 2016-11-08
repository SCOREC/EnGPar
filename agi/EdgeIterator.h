#ifndef __EDGE_ITERATOR__H__
#define __EDGE_ITERATOR__H__
#include "agi.h"

namespace agi {

class Ngraph;

class EdgeIterator {
  friend class Ngraph;
 private:
  etype type;
  lid_t* loc;
  EdgeIterator(etype t, lid_t* l) : type(t), loc(l) {}
};
  
}


#endif
