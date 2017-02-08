#ifndef ENGPAR_VTXBOUNDS_H
#define ENGPAR_VTXBOUNDS_H

namespace engpar {
  //class VtxBounds : public engpar::Bounds
  class VtxBounds  {
  public:
    VtxBounds(agi::Ngraph* g) {
      total_boundaries =0;
      init(g);
    }
  private:
    void init(agi:Ngraph* g) {
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g->begin(0);
      while ((edge = g->iterate(eitr))) {
	agi::GraphVertex* pin;
	agi::PinIterator* pitr = g->pins(edge);
	int deg = g->degree(edge);
	for (int i=0;i<deg;i++) {
	  pin = g->iterate(pitr);
	  if (PCU_Comm_Self()!=g->owner(pin)){ 
	  }
	}
      }
    }
  protected:
    std::unordered_map<int,int> b;
    int total_boundaries;
    
  }
}
//engpar::Bounds makeVtxBounds(agi::Ngraph graph) {
engpar::VtxBounds makeVtxBounds(agi::Ngraph graph) {
  return new VtxBounds(graph);
}

#endif
