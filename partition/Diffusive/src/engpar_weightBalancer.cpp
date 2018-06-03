#include "engpar_balancer.h"
#include "engpar_weightSelector.h"
#include "engpar_queue.h"
#include <PCU.h>
#include "../engpar.h"
#include <engpar_support.h>

namespace engpar {


  class WeightBalancer : public Balancer {
  public:
    WeightBalancer(Input* input_, int v);
    ~WeightBalancer() {}

    bool runStep(double tol);
    void balance();
  };
  
  WeightBalancer::WeightBalancer(Input* input_, int v) :
    Balancer(input_,v,"weightBalancer") {
  }
  
  bool WeightBalancer::runStep(double tolerance) {
    WeightInput* inp = dynamic_cast<WeightInput*>(input);
    double stepTime = PCU_Time();
    double imb = EnGPar_Get_Imbalance(getWeight(input->g,-1,false));
    //Check for completition of criteria
    if (imb < tolerance)
      return false;
    //Check stagnation detection
    sd->push(imb);
    if (sd->isFull()&&sd->slope()>0)
      return false;
    Sides* sides = makeSides(inp->g,inp->primary_edge_type);
    if (verbosity>=3)
      EnGPar_Status_Message("%d: %s\n",PCU_Comm_Self(), sides->print("Sides").c_str());
    Weights* targetWeights = makeWeights(inp->g,false, sides,-1);
    if (verbosity>=3)
      EnGPar_Status_Message("%d: %s\n",PCU_Comm_Self(),
             targetWeights->print("Weights").c_str());
    Weights** completedWs = NULL;
    Targets* targets = makeTargets(inp->g->isHyper(),inp->step_factor,sides,targetWeights,
                                   sideTol, completedWs,completed_weights);
    delete sides;
    if (completedWs) {
      for (unsigned int i=0;i<completed_dimensions.size();i++)
        delete completedWs[i];
      delete [] completedWs;
    }
    delete targetWeights;
        
    if (verbosity>=3)
      EnGPar_Status_Message("%d: %s\n",PCU_Comm_Self(), targets->print("Targets").c_str());
    Queue* pq;
    double t = PCU_Time();
    if (false)
      //TODO: Come up with potential iteration queues for weight diffusion
      //pq = createDistanceQueue(inp);
      ;
    else 
      pq = createIterationQueue(input->g);
    distance_time+=PCU_Time()-t;

    WeightSelector* selector = makeWeightSelector(inp,pq);
    agi::WeightMigration* plan = new agi::WeightMigration(input->g);
    wgt_t planW = 0.0;
    int selectIterations=10;
    for (int i=0;i<selectIterations;i++) {
      wgt_t beforeW = planW;
      planW = selector->select(targets,plan,beforeW);
      if (planW-beforeW<.00001) {
        break;
      }
    }
    delete pq;
    delete targets;
    delete selector;
    
    stepTime = PCU_Time()-stepTime;
    int numMigrate = plan->size();
    numMigrate = PCU_Add_Int(numMigrate);

    if (numMigrate>0)
      input->g->migrate(plan, migrTime);
    else {
      //delete plan;
    }
    if (verbosity >= 1) {
      if (!PCU_Comm_Self()) {
        EnGPar_Status_Message("  Step took %f seconds\n",stepTime);
        EnGPar_Status_Message("  Imbalances <v, e0, ...>: ");
      }
      printImbalances(input->g);
      totStepTime+=stepTime;
    }
    if (verbosity >= 2) {
      if (!PCU_Comm_Self()) {
        if (sd->isFull())
          EnGPar_Status_Message("    Slope: %f\n",sd->slope());
        EnGPar_Status_Message("    Migrating %d weight took %f seconds\n",numMigrate, 0.0);//migrTime->getTime("total"));
      }
    }

    if (numMigrate == 0)
      return false;

    return true; //not done balancing
  }
  void WeightBalancer::balance() {
    WeightInput* inp = dynamic_cast<WeightInput*>(input);
    if (EnGPar_Is_Log_Open()) {
      /*
        TODO: write logging info for weight input
      char message[1000];
      sprintf(message,"balance() : \n");
      //Log the input parameters
      sprintf(message,"%s priorities :",message);
      for (unsigned int i=0;i<inp->priorities.size();i++)
        sprintf(message,"%s %d",message,inp->priorities[i]);
      sprintf(message,"%s\n",message);
      if (inp->tolerances.size()>0) {
        sprintf(message,"%s tolerances :",message);
        for (unsigned int i=0;i<inp->tolerances.size();i++)
          sprintf(message,"%s %f",message,inp->tolerances[i]);
        sprintf(message,"%s\n",message);
      }
      sprintf(message,"%s maxIterations : %d\n",message,inp->maxIterations);
      sprintf(message,"%s maxIterationsPerType : %d\n",
              message,inp->maxIterationsPerType);
      sprintf(message,"%s step_factor : %f\n",message,inp->step_factor);
      sprintf(message,"%s sides_edge_type : %d\n",message,inp->sides_edge_type);
      sprintf(message,"%s selection_edge_type : %d\n",
              message,inp->selection_edge_type);
      sprintf(message,"%s countGhosts : %d\n",message,inp->countGhosts);
      
      EnGPar_Log_Function(message);
      */
    }
    if (1 == PCU_Comm_Peers()) {
      EnGPar_Warning_Message("Ran in serial, nothing to do exiting...\n");
      return;
    }

    //TODO:check to make sure these timers still make sense
    //Setup the migration timers
    // migrTime->addTimer("setup");
    // migrTime->addTimer("comm");
    // migrTime->addTimer("build");
    // migrTime->addTimer("total");

    //TODO:replace original owners with original weights
    //Setup the original owners arrays before balancing
    // input->g->setOriginalOwners();
    // unsigned int index=0;
    // target_dimension = inp->priorities[index];
    
    //Construct Stagnation detection
    sd = new SDSlope;
    

    //Set side tolerance to arbitrarily high number such that it will be ignored in targeting
    sideTol = inp->g->numGlobalEdges()*10;
    
    if (!PCU_Comm_Self() && verbosity >= 0)
      EnGPar_Status_Message("Starting weight diffusion with imbalance: ");
    printImbalances(input->g);

    int step = 0;
    double time = PCU_Time();
    //double targetTime=PCU_Time();
    while (step++<inp->maxIterations&&runStep(inp->tol));
    delete sd;
    time = PCU_Time()-time;

    if (verbosity >= 0) {
      time = PCU_Max_Double(time);
      if (!PCU_Comm_Self()) {
        if(step==inp->maxIterations)
          EnGPar_Status_Message("EnGPar ran to completion in %d iterations in %f seconds\n",
                 inp->maxIterations, time);
        else
          EnGPar_Status_Message("EnGPar converged in %d iterations in %f seconds\n",
                 step,time);
      }
    }
    /*
    if (verbosity >= 2) {
      printMigrationStats(migrTime);
      double maxMigr = migrTime->processMax("total");
      double maxPlan = PCU_Max_Double(totStepTime);
      distance_time = PCU_Max_Double(distance_time);
      if (!PCU_Comm_Self()) {
        EnGPar_Status_Message("Migration took %f s, %f%% of the total time\n", maxMigr, maxMigr/time*100);
        EnGPar_Status_Message("Planning took %f s, %f%% of the total time\n", maxPlan, maxPlan/time*100);
        EnGPar_Status_Message("Distance Computation (part of Planning) took %f seconds, %f%% of the total time\n",
            distance_time, distance_time/time*100);
      }
    }
    */
    if (EnGPar_Is_Log_Open())
      EnGPar_End_Function();
  }

  void balanceWeights(WeightInput* in, int v_) {
    if (EnGPar_Is_Log_Open()) {
      char message[25];
      sprintf(message,"balanceWeights\n");
      EnGPar_Log_Function(message);
      EnGPar_End_Function();
    }
    WeightBalancer* balancer = new WeightBalancer(in,v_);
    balancer->balance();
    delete balancer;
  }
}

