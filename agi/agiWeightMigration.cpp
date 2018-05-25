#include "agiMigration.h"
#include "ngraph.h"
#include <PCU.h>
namespace agi {
  WeightMigration::WeightMigration(Ngraph* graph) {
    g = graph;
    s = 0;
  }
  WeightMigration::~WeightMigration() {
  }

  lid_t WeightMigration::size() const {
    return s;
  }
  bool WeightMigration::has(GraphVertex* v) const {
    return sending.find(v)!=sending.end();
  }
  bool WeightMigration::has(GraphVertex* v, GraphVertex* u) const {
    WMPlan::const_iterator itr = sending.find(v);
    return itr!=sending.end()&& itr->second.find(u) != itr->second.end();
  }
  wgt_t WeightMigration::get(GraphVertex* v,GraphVertex* u) const {
    return sending.find(v)->second.find(u)->second;
  }

  void WeightMigration::insert(GraphVertex* v, GraphVertex* u, wgt_t w) {
    sending[v][u]+=w;
    s+=w;
  }

  void WeightMigration::clear() {
    sending.clear();
  }

  WeightMigration::iterator WeightMigration::begin() const {
    return iterator(&sending,sending.begin());
  }
  WeightMigration::iterator WeightMigration::end() const {
    return iterator(&sending,sending.end());
  }
}
