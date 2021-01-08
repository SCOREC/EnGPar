#ifndef ENGPAR_BALANCER_H
#define ENGPAR_BALANCER_H

#include "engpar_sides.h"
#include "engpar_weights.h"
#include "engpar_targets.h"
#include "engpar_selector.h"
#include "engpar_sd.h"
#include "agiMigrationTimers.h"
#include "../engpar_diffusive_input.h"
#include <torch/script.h>

namespace engpar {
  wgt_t getMaxWeight(agi::Ngraph*, int);
  wgt_t getAvgWeight(agi::Ngraph*, int);
  class Balancer {
  public:
    Balancer(agi::Ngraph*& graph_, double factor_, int verbosity_,
             const char* name_);
    Balancer(Input* input_,int verbosity_,const char* name_);
    virtual ~Balancer();
    virtual int runStep(double tolerance);
    int runStepTorch(double tolerance, bool useTorch, torch::jit::script::Module& module);
    void balance();
    void partWeightBalancer(Sides* sides, double tolerance);
    
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
    agi::Ngraph* weightGraph;
    agi::MigrationTimers* migrTime;
  };

  class WeightBalancer : public Balancer {
  public:
    WeightBalancer(Input* input_, int v);
    ~WeightBalancer() {}

    int runStep(double tol);
    void balance();
  };

}

#endif
