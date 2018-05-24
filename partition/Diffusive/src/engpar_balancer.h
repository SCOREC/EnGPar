#ifndef ENGPAR_BALANCER_H
#define ENGPAR_BALANCER_H

#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"
#include "engpar_selector.h"
#include "engpar_sd.h"
#include "agiMigrationTimers.h"

namespace engpar {
  wgt_t getMaxWeight(agi::Ngraph*, int);
  wgt_t getAvgWeight(agi::Ngraph*, int);
  class Balancer {
  public:
    Balancer(agi::Ngraph*& graph_, double factor_, int verbosity_,
             const char* name_);
    Balancer(Input* input_,int verbosity_,const char* name_);
    virtual ~Balancer() {delete input; delete migrTime;}
    virtual bool runStep(double tolerance);
    void balance();
    
  protected:
    agi::Ngraph* graph;
    int verbosity;
    const char* name;

    int target_dimension;
    SDSlope* sd;
    std::vector<int> completed_dimensions;
    std::vector<wgt_t> completed_weights;
    Input* input;
    double totStepTime;
    double distance_time;
    int sideTol;
    agi::MigrationTimers* migrTime;
  };
}

#endif
