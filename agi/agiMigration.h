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

    void insert(std::pair<GraphVertex*, int>);

    void clear();
    typedef std::vector<GraphVertex*>::iterator iterator;
    iterator begin();
    iterator end();
  private:
    Ngraph* g;
    GraphTag* dest;
    std::vector<GraphVertex*> sending;
  };
  
}
#endif
