#include "agiMigration.h"
#include "ngraph.h"
#include <PCU.h>
namespace agi {
  Migration::Migration(Ngraph* graph) {
    g = graph;
    dest = NULL;
    if (g->numLocalVtxs()!=0) {
      dest = g->createLongTag();
      VertexIterator* vitr = g->begin();
      GraphVertex* v;
      while ((v = g->iterate(vitr))) {
        g->setLongTag(dest,v,-1);
      }
    }
  }
  Migration::~Migration() {
    if (dest)
      g->destroyTag(dest);
  }

  lid_t Migration::size() const {
    return sending.size();
  }
  bool Migration::has(GraphVertex* v) const {
    return g->localID(v)<g->numLocalVtxs()&&g->getLongTag(dest,v) !=-1;
  }
  int Migration::get(GraphVertex* v) const {
    return g->getLongTag(dest,v);
  }

  void Migration::insert(std::pair<GraphVertex*, int> pair) {
    if (has(pair.first)) {
      g->setLongTag(dest,pair.first,pair.second);
      return;
    }
    sending.push_back(pair.first);
    g->setLongTag(dest,pair.first,pair.second);
  }

  void Migration::clear() {
    for (unsigned int i=0;i<sending.size();i++) {
      g->setLongTag(dest,sending[i],-1);
    }
    sending.clear();
  }

  Migration::iterator Migration::begin() {
    return sending.begin();
  }
  Migration::iterator Migration::end() {
    return sending.end();
  }
}
