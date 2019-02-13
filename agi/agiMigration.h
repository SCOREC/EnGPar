#ifndef __AGI__MIGRATION_H__
#define __AGI__MIGRATION_H__
#include "agi.h"
#include <vector>

namespace agi {

  class Ngraph;
  class GraphTag;

  class Migration {
  public:
    Migration(Ngraph*);
    ~Migration();

    lid_t size() const;

    bool has(GraphVertex*) const;
    int get(GraphVertex*) const;

    bool insert(std::pair<GraphVertex*, int>);

    void clear();
    typedef std::vector<GraphVertex*>::iterator iterator;
    iterator begin();
    iterator end();
  private:
    Ngraph* g;
    GraphTag* dest;
    std::vector<GraphVertex*> sending;
  };

  //Map from local vtx -> ghost vtx -> weight to send
  typedef std::unordered_map<agi::GraphVertex*, wgt_t> WMData;
  typedef std::unordered_map<agi::GraphVertex*, WMData > WMPlan;

  
  class WMIterator {
  public:
    WMIterator(const WMPlan* p,WMPlan::const_iterator start) {
      plan = p;
      outer_itr = start;
      if (outer_itr!=plan->end())
        inner_itr = outer_itr->second.begin();
    }
    GraphVertex* u() {
      return outer_itr->first;
    }
    GraphVertex* v() {
      return inner_itr->first;
    }
    wgt_t weight() {
      return inner_itr->second;
    }
    WMIterator& operator++() {
      inner_itr++;
      if (inner_itr==outer_itr->second.end()) {
        outer_itr++;
        if (outer_itr!=plan->end())
          inner_itr = outer_itr->second.begin();
      }
      return *this;
    }
    WMIterator operator++(int) {
      WMIterator copy = *this;
      inner_itr++;
      if (inner_itr==outer_itr->second.end()) {
        outer_itr++;
        if (outer_itr!=plan->end())
          inner_itr = outer_itr->second.begin();
      }
      return copy;
      
    }
    bool operator==(const WMIterator& other) {
      return outer_itr == other.outer_itr &&
            (outer_itr==plan->end() || inner_itr == other.inner_itr);
    }
    bool operator!=(const WMIterator& other) {
      return !(*this == other);
    }
  private:
    const WMPlan* plan;
    WMPlan::const_iterator outer_itr;
    WMData::const_iterator inner_itr;
    
  };

  class WeightMigration {
  public:
    WeightMigration(Ngraph*);
    ~WeightMigration();
    lid_t size() const;

    bool  has(GraphVertex*) const;
    bool  has(GraphVertex*, GraphVertex*) const;
    wgt_t get(GraphVertex*, GraphVertex*) const;

    void insert(GraphVertex* v, GraphVertex* u, wgt_t w);

    void clear();
    typedef WMIterator iterator;
    iterator begin() const;
    iterator end() const;

  private:
    Ngraph* g;
    wgt_t s;
    WMPlan sending;
  };
}
#endif
