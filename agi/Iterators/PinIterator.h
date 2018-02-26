#ifndef __PIN_ITERATOR__H__
#define __PIN_ITERATOR__H__
#include "../agi.h"
#include <stdint.h>
namespace agi {

class Ngraph;

//A class to handle iteration over edges of a certain type
//The pointer = 1+type+num_types*edge_id
class PinIterator {
  friend class Ngraph;
  lid_t* loc;
  lid_t* end;
  PinIterator(lid_t* l,lid_t* e)
    :loc(l),end(e) {}
  void iterate() {loc++;}
};
  
}


#endif
