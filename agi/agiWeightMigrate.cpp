#include "ngraph.h"
#include <PCU.h>
#include <unordered_set>
#include <vector>
#include <engpar_support.h>
#include "agiMigrationTimers.h"


namespace agi {

  void Ngraph::migrate(WeightMigration* plan, MigrationTimers* mt) {
    delete plan;
  }

}
