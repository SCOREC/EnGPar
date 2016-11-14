#ifndef EDGE_H__
#define EDGE_H__
#include "agi.h"

namespace agi {

class Ngraph;
class apfGraph;
  
class Edge {
  friend class Ngraph;
  friend class apfGraph;
 private:
  lid_t lid;
  gid_t gid;
  wgt_t weight;
  etype type;
  Edge() {}
  Edge(lid_t l,gid_t g,wgt_t w,etype t) : lid(l), gid(g), weight(w), type(t) {};
};

}

#endif
