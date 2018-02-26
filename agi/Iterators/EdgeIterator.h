#ifndef __EDGE_ITERATOR__H__
#define __EDGE_ITERATOR__H__
#include "../agi.h"
#include <stdint.h>
namespace agi {

class Ngraph;

//A class to handle iteration over edges of a certain type
//The pointer = 1+type+num_types*edge_id
class EdgeIterator {
  friend class Ngraph;
 public:
  virtual ~EdgeIterator() {}
 protected:
  lid_t* loc;
  lid_t* end;
  int num_types;
  EdgeIterator(etype t, int nt, lid_t* l,lid_t deg)
    :loc((lid_t*)(1+t+nt*(uintptr_t)l)),
    end((lid_t*)(1+t+nt*((uintptr_t)l+deg))),
    num_types(nt) {}
  void iterate() {loc = (lid_t*)((uintptr_t)loc+num_types);}
  virtual bool isHyper() {return false;}
};
  
}


#endif
