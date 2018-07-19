#include "engpar_balancer.h"
#include "engpar_queue.h"
#include <PCU.h>
#include <engpar_metrics.h>
#include <engpar_support.h>
#include <agiMigration.h>
namespace {
  void printMigrationStats(agi::MigrationTimers* migrTime) {
    double maxSetup = migrTime->processMax("setup");
    double maxComm = migrTime->processMax("comm");
    double maxBuild = migrTime->processMax("build");
    double maxTot = migrTime->processMax("total");
    double localTime[4];
    localTime[0] = migrTime->getTime("setup");
    localTime[1] = migrTime->getTime("comm");
    localTime[2] = migrTime->getTime("build");
    localTime[3] = migrTime->getTime("total");
    double ratios[4], globalRatios[4];
    ratios[0] = localTime[0]/localTime[3]; // setup/total
    ratios[1] = localTime[1]/localTime[3]; // comm/total
    ratios[2] = localTime[2]/localTime[3]; // build/total
    ratios[3] = (localTime[0]+localTime[1]+localTime[2])/localTime[3]; // (setup+comm+build)/total
    for(int i=0; i<4; i++) globalRatios[i] = ratios[i];
    PCU_Max_Doubles(globalRatios,4);
    int count = migrTime->getCount("comm");
    if (!PCU_Comm_Self() && count) {
      EnGPar_Status_Message("max migration time (s) "
          "<total, setup, comm, build> = <%f, %f, %f, %f>\n",
          maxTot, maxSetup, maxComm, maxBuild);
      EnGPar_Status_Message("max migration ratios "
          "<setup/total, comm/total, build/total, (setup+comm+build)/total> = <%f, %f, %f, %f>\n",
          globalRatios[0], globalRatios[1], globalRatios[2], globalRatios[3]);
    }
    for(int i=0; i<4; i++) globalRatios[i] = ratios[i];
    PCU_Min_Doubles(globalRatios,4);
    if (!PCU_Comm_Self() && count) {
      EnGPar_Status_Message("min migration ratios "
          "<setup/total, comm/total, build/total, (setup+comm+build)/total> = <%f, %f, %f, %f>\n",
          globalRatios[0], globalRatios[1], globalRatios[2], globalRatios[3]);
    }
  }
}

namespace engpar {

  wgt_t getMaxWeight(agi::Ngraph* g, int dimension, bool countGhosts) {
    wgt_t w = getWeight(g,dimension,countGhosts);
    return PCU_Max_Double(w);
  }

  wgt_t getAvgWeight(agi::Ngraph* g, int dimension, bool countGhosts) {
    wgt_t w = getWeight(g,dimension,countGhosts);
    return PCU_Add_Double(w) / PCU_Comm_Peers();
  }
  double averageSides(Sides* s) {
    double tot = s->total();
    tot = PCU_Add_Double(tot);
    return tot / PCU_Comm_Peers();
  }
  
  Balancer::Balancer(agi::Ngraph*& g, double f, int v, const char* n) :
    graph(g), verbosity(v), name(n) {
    input = createDiffusiveInput(g,f);
    totStepTime=0;
    distance_time=0;
    migrTime = new agi::MigrationTimers;
    weightGraph = NULL;
  }
  Balancer::Balancer(Input* input_, int v, const char* n) :
    graph(input_->g),verbosity(v), name(n), input(input_) {
    totStepTime=0;
    distance_time=0;
    migrTime = new agi::MigrationTimers;
    weightGraph = NULL;
  }
  Balancer::~Balancer() {
    delete input;
    delete migrTime;
    if (weightGraph)
      agi::destroyGraph(weightGraph);
  }

  bool Balancer::runStep(double tolerance) {
    DiffusiveInput* inp = dynamic_cast<DiffusiveInput*>(input);
    double stepTime = PCU_Time();
    double imb = EnGPar_Get_Imbalance(getWeight(input->g,target_dimension,inp->countGhosts));
    //Check for completition of criteria
    if (imb < tolerance)
      return false;
    //Check stagnation detection
    sd->push(imb);
    if (sd->isFull()&&sd->slope()>0)
      return false;
    Sides* sides = makeSides(inp);
    if (verbosity>=3)
      EnGPar_Status_Message("%d: %s\n",PCU_Comm_Self(), sides->print("Sides").c_str());
    Weights* targetWeights = makeWeights(inp, sides,target_dimension);
    if (verbosity>=3)
      EnGPar_Status_Message("%d: %s\n",PCU_Comm_Self(),
             targetWeights->print("Weights").c_str());
    Weights** completedWs = NULL;
    if (completed_dimensions.size()>0) {
      completedWs= new Weights*[completed_dimensions.size()];
      for (unsigned int i=0;i<completed_dimensions.size();i++) {
        completedWs[i] = makeWeights(inp,sides,completed_dimensions[i]);
      }
    }
    Targets* targets;
    if (inp->runPartWeightBalancer) {
      targets = makePartWeightTargets(inp, sides, weightGraph->getWeightPartition(), sideTol,
                                      completedWs, completed_weights);
    }
    else {
      targets = makeTargets(inp,sides,targetWeights,sideTol,
                                     completedWs,completed_weights);
    }
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
    if (inp->useDistanceQueue) {
      pq = createDistanceQueue(inp);
    }
    else 
      pq = createIterationQueue(input->g);
    distance_time+=PCU_Time()-t;
    Selector* selector = makeSelector(inp,pq,&completed_dimensions,
                                      &completed_weights);
    agi::Migration* plan = new agi::Migration(input->g);
    wgt_t planW = 0.0;
    for (unsigned int cavSize=2;cavSize<=12;cavSize+=2) {
      planW = selector->select(targets,plan,planW,cavSize,target_dimension);
    }
    selector->selectDisconnected(plan,target_dimension);
    if (completed_dimensions.size()>0) {
      int sizes[2];
      sizes[0] = plan->size();
      Midd* midd = selector->trim(targets,plan);
      selector->cancel(plan,midd);
      sizes[1]=plan->size();
      if (verbosity>=2) {
        PCU_Add_Ints(sizes,2);
        if (!PCU_Comm_Self())
          EnGPar_Status_Message("  Plan was trimmed from %d to %d vertices\n",sizes[0],sizes[1]);
      }
    }
    delete pq;
    delete targets;

    if (inp->runPartWeightBalancer) {
      selector->updatePartWeight(plan,target_dimension, weightGraph->getWeightPartition());
    }
    delete selector;

    stepTime = PCU_Time()-stepTime;
    int numMigrate = plan->size();
    numMigrate = PCU_Add_Int(numMigrate);

    if (numMigrate>0)
      input->g->migrate(plan, migrTime);
    else
      delete plan;
    
    if (verbosity >= 1) {
      char buffer[100];
      getImbalances(input->g,buffer);
      if (!PCU_Comm_Self()) {
        EnGPar_Status_Message("  Step took %f seconds\n",stepTime);
        EnGPar_Status_Message("  Imbalances <v, e0, ...>: %s\n",buffer);
      }
      totStepTime+=stepTime;
    }
    if (verbosity >= 2) {
      if (!PCU_Comm_Self()) {
        if (sd->isFull())
          EnGPar_Status_Message("    Slope: %f\n",sd->slope());
        EnGPar_Status_Message("    Migrating %d vertices took %f seconds\n",numMigrate, migrTime->getTime("total"));
      }
    }

    if (numMigrate == 0)
      return false;

    return true; //not done balancing
  }
  void Balancer::balance() {
    DiffusiveInput* inp = dynamic_cast<DiffusiveInput*>(input);
    if (EnGPar_Is_Log_Open()) {
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
    }
    if (1 == PCU_Comm_Peers()) {
      EnGPar_Warning_Message("ran in serial, nothing to do exiting...\n");
      return;
    }

    //Setup the migration timers
    migrTime->addTimer("setup");
    migrTime->addTimer("comm");
    migrTime->addTimer("build");
    migrTime->addTimer("total");

    //Setup the original owners arrays before balancing
    input->g->setOriginalOwners();
    unsigned int index=0;
    target_dimension = inp->priorities[index];
    
    //Construct Stagnation detection
    sd = new SDSlope;
    
    //Set imbalance tolerance
    double tol=1.1;
    if (inp->tolerances.size()>index)
      tol = inp->tolerances[index];
    
    //Set side tolerance
    Sides* sides = makeSides(inp);
    sideTol = averageSides(sides);

    if (inp->runPartWeightBalancer) {
      if (!PCU_Comm_Self())
        EnGPar_Status_Message("Starting part weight balancer on type %d\n",target_dimension);
      partWeightBalancer(sides,tol);
    }
    delete sides;
    
    if (verbosity >= 0) {
      char buffer[100];
      getImbalances(input->g, buffer);
      if (!PCU_Comm_Self())
        EnGPar_Status_Message("Starting criteria type %d with imbalances: %s\n",
                              target_dimension,buffer);
    }
    if (!PCU_Comm_Self() && verbosity >= 1)
      EnGPar_Status_Message("Side Tolerance is: %d\n", sideTol);

    int step = 0;
    int inner_steps=0;
    double time = PCU_Time();
    double targetTime=PCU_Time();
    while (step++<inp->maxIterations) {
      //runStep(tol) balances the current dimension.
      //Advance to the next dimension if the current dimesion is balanced
      //or the maximum per dimension iterations is reached.
      if (!runStep(tol) || inner_steps++ >= inp->maxIterationsPerType) {
        // Set the imbalance limit for the higher priority dimension
        // while balancing lower priority dimensions to be the larger of
        // the specified imbalance (tgtMaxW) or, if it wasn't reached, the
        // current weight (maxW).
        completed_dimensions.push_back(target_dimension);
        double maxW = getMaxWeight(input->g,target_dimension, inp->countGhosts);
        double tgtMaxW = getAvgWeight(input->g,target_dimension,inp->countGhosts)*tol;
        maxW = ( maxW < tgtMaxW ) ? tgtMaxW : maxW;
        completed_weights.push_back(maxW);
        
        targetTime = PCU_Time()-targetTime;
        targetTime = PCU_Max_Double(targetTime);
        if (verbosity >= 0 && !PCU_Comm_Self()) {
          EnGPar_Status_Message("Completed criteria type %d in %d steps and took %f seconds\n",
                 target_dimension, inner_steps, targetTime);
        }
        targetTime=PCU_Time();
        
        index++;
        if (index==inp->priorities.size())
          break;
        inner_steps=0;
        //Set new target criteria
        target_dimension=inp->priorities[index];
        //Recreate stagnation detection
        delete sd;
        sd = new SDSlope;
        //Set new tolerance
        if (inp->tolerances.size()>index)
          tol = inp->tolerances[index];        
        //Set side tolerance
        Sides* sides = makeSides(inp);
        sideTol = averageSides(sides);

        if (inp->runPartWeightBalancer) {
          if (!PCU_Comm_Self())
            EnGPar_Status_Message("Starting part weight balancer on type %d\n",target_dimension);
          partWeightBalancer(sides,tol);
        }

        delete sides;

        char buffer[100];
        getImbalances(input->g,buffer);
        if (!PCU_Comm_Self() && verbosity >= 0)
          EnGPar_Status_Message("Starting criteria type %d with imbalances: %s\n",
                                target_dimension,buffer);

        if (!PCU_Comm_Self() && verbosity >= 1)
          EnGPar_Status_Message("Side Tolerance is: %d\n", sideTol);


      }      
    }
    delete sd;
    time = PCU_Time()-time;


    if (verbosity >= 0) {
      time = PCU_Max_Double(time);
      if (!PCU_Comm_Self()) {
        if(step==inp->maxIterations)
          EnGPar_Status_Message("EnGPar ran to completion in %d iterations in %f seconds\n",
                 inp->maxIterations, time);
        else
          EnGPar_Status_Message("EnGPar converged in %lu iterations in %f seconds\n",
                 step-inp->priorities.size(),time);
      }
    }
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
    if (EnGPar_Is_Log_Open())
      EnGPar_End_Function();
  }
  void Balancer::partWeightBalancer(Sides* sides,double tol) {
    DiffusiveInput* inp = dynamic_cast<DiffusiveInput*>(input);

    if (weightGraph)
      agi::destroyGraph(weightGraph);
    weightGraph = agi::createEmptyGraph();
    agi::lid_t num_verts = 1;
    agi::gid_t* verts = new agi::gid_t[num_verts];
    agi::wgt_t* weights = new agi::wgt_t[num_verts];
    verts[0] = PCU_Comm_Self();
    weights[0] = getWeight(input->g,target_dimension,inp->countGhosts);
    agi::lid_t num_edges = sides->size();
    agi::lid_t* degrees = new agi::lid_t[num_edges];
    agi::gid_t* edges = new agi::gid_t[num_edges];
    agi::gid_t* pins = new agi::gid_t[num_edges*2];
    agi::lid_t num_ghosts = num_edges;
    agi::gid_t* ghosts = new agi::gid_t[num_edges];
    agi::part_t* owners = new agi::part_t[num_edges];
    Sides::iterator itr;
    int i = 0;
    for (itr = sides->begin(); itr != sides->end(); itr++, ++i) {
      degrees[i] = 2;
      edges[i] = i;
      pins[2*i] = PCU_Comm_Self();
      pins[2*i+1] = itr->first;
      ghosts[i] = itr->first;
      owners[i] = itr->first;
    }
    weightGraph->constructVerts(false,num_verts,verts,weights);
    weightGraph->constructEdges(num_edges,edges,degrees,pins);
    weightGraph->constructGhosts(num_ghosts,ghosts,owners);

    //Set the tolerance of the part weight balancer slightly less to give some wiggle room
    //TODO: experiment on this `.75` parameter
    WeightInput* input = createWeightInput(weightGraph,tol*.75,inp->step_factor);

    WeightBalancer* weightBalancer = new WeightBalancer(reinterpret_cast<Input*>(input),-1);
    weightBalancer->balance();
    delete weightBalancer;
    double imb = EnGPar_Get_Imbalance(getWeight(weightGraph,-1,false));
    if (!PCU_Comm_Self()) {
      EnGPar_Status_Message("Part Weight Balancer finished with imbalance %.3f\n",imb);
    }

    delete [] verts;
    delete [] weights;
    delete [] edges;
    delete [] degrees;
    delete [] pins;
    delete [] ghosts;
    delete [] owners;
  }

  void balance(Input* in, int v_) {
    if (EnGPar_Is_Log_Open()) {
      char message[25];
      sprintf(message,"makeBalancer\n");
      EnGPar_Log_Function(message);
      EnGPar_End_Function();
    }
    Balancer* balancer = new Balancer(in,v_,"balancer");
    balancer->balance();
    delete balancer;
  }
}

