#include "engpar_balancer.h"
#include "engpar_queue.h"
#include <PCU.h>
#include "../engpar.h"
#include <engpar_support.h>
namespace engpar {

  wgt_t getMaxWeight(agi::Ngraph* g, int dimension) {
    wgt_t w = getWeight(g,dimension);
    return PCU_Max_Double(w);
  }

  wgt_t getAvgWeight(agi::Ngraph* g, int dimension) {
    wgt_t w = getWeight(g,dimension);
    return PCU_Add_Double(w) / PCU_Comm_Peers();
  }
  
  Balancer::Balancer(agi::Ngraph* g, double f, int v, const char* n) :
    agi::Balancer(g,v,n) {
    input = new Input(g);
    input->step_factor = f;
    times[0]=0;
    times[1]=0;
    distance_time=0;
  }
  Balancer::Balancer(Input* input_, int v, const char* n) :
    agi::Balancer(input_->g,v,n), input(input_) {
    times[0]=0;
    times[1]=0;
    distance_time=0;
  }
  bool Balancer::runStep(double tolerance) {
    double time[2];
    time[0] = PCU_Time();

    double imb = EnGPar_Get_Imbalance(getWeight(input->g,target_dimension));
    //Check for completition of criteria
    if (imb < tolerance)
      return false;
    //Check stagnation detection
    sd->push(imb);
    if (sd->isFull()&&sd->slope()>0)
      return false;
    Sides* sides = makeSides(input);
    if (verbosity>=3)
      printf("%d: %s\n",PCU_Comm_Self(), sides->print("Sides").c_str());
    Weights* targetWeights = makeWeights(input, sides,target_dimension);
    if (verbosity>=3)
      printf("%d: %s\n",PCU_Comm_Self(),
             targetWeights->print("Weights").c_str());
    Weights** completedWs = NULL;
    if (completed_dimensions.size()>0) {
      completedWs= new Weights*[completed_dimensions.size()];
      for (unsigned int i=0;i<completed_dimensions.size();i++) {
        completedWs[i] = makeWeights(input,sides,completed_dimensions[i]);
      }
    }
    Targets* targets = makeTargets(input,sides,targetWeights,
                                   completedWs,completed_weights);
    delete sides;
    if (completedWs) {
      for (unsigned int i=0;i<completed_dimensions.size();i++)
        delete completedWs[i];
      delete [] completedWs;
    }
    delete targetWeights;
        
    if (verbosity>=3)
      printf("%d: %s\n",PCU_Comm_Self(), targets->print("Targets").c_str());
    Queue* pq;
    double t = PCU_Time();
    if (input->useDistanceQueue) {
      pq = createDistanceQueue(input->g);
    }
    else 
      pq = createIterationQueue(input->g);
    distance_time+=PCU_Time()-t;
    Selector* selector = makeSelector(input,pq,&completed_dimensions,
                                      &completed_weights);
    agi::Migration* plan = new agi::Migration(input->g);
    wgt_t planW = 0.0;
    for (unsigned int cavSize=2;cavSize<=12;cavSize+=2) {
      planW += selector->select(targets,plan,planW,cavSize,target_dimension);
    }
    if (completed_dimensions.size()>0) {
      int sizes[2];
      sizes[0] = plan->size();
      Midd* midd = selector->trim(targets,plan);
      selector->cancel(plan,midd);
      sizes[1]=plan->size();
      if (verbosity>=2) {
        PCU_Add_Ints(sizes,2);
        if (!PCU_Comm_Self())
          printf("Plan was trimmed from %d to %d vertices\n",sizes[0],sizes[1]);
      }
    }
    delete pq;
    delete targets;
    delete selector;
    
    time[0] = PCU_Time()-time[0];
    int numMigrate = plan->size();
    numMigrate = PCU_Add_Int(numMigrate);
    if (verbosity>=3) {
      int* counts = new int[PCU_Comm_Peers()];
      for (int i=0;i<PCU_Comm_Peers();i++)
        counts[i] = 0;
      agi::Migration::iterator itr;
      for (itr = plan->begin();itr!=plan->end();itr++)
        counts[plan->get(*itr)]++;
      for (int i=0;i<PCU_Comm_Peers();i++)
        if (counts[i]>0)
          printf("%d sending %d to %d\n",PCU_Comm_Self(),counts[i],i);
      delete [] counts;
    }

    time[1] = PCU_Time();
    if (numMigrate>0)
      input->g->migrate(plan);
    else
      delete plan;
    time[1] = PCU_Time()-time[1];
    
    if (verbosity >= 1) {
      if (!PCU_Comm_Self()) {
        printf("Step took %f seconds\n",time[0]);
        printf("Imbalances <v, e0, ...>: ");
      }
      printImbalances(input->g);
      times[0]+=time[0];      
    }
    if (verbosity >= 2) {
      if (!PCU_Comm_Self()) {
        if (sd->isFull())
          printf("Slope: %f\n",sd->slope());
        printf("Migrating %d vertices took %f seconds\n",numMigrate,time[1]);
      }
      times[1]+=time[1];
    }

    if (numMigrate == 0)
      return false;

    return true; //not done balancing
  }
  void Balancer::balance(double) {
    if (EnGPar_Is_Log_Open()) {
      char message[1000];
      sprintf(message,"balance() : \n");
      //Log the input parameters
      sprintf(message,"%s priorities :",message);
      for (unsigned int i=0;i<input->priorities.size();i++)
        sprintf(message,"%s %d",message,input->priorities[i]);
      sprintf(message,"%s\n",message);
      if (input->tolerances.size()>0) {
        sprintf(message,"%s tolerances :",message);
        for (unsigned int i=0;i<input->tolerances.size();i++)
          sprintf(message,"%s %f",message,input->tolerances[i]);
        sprintf(message,"%s\n",message);
      }
      sprintf(message,"%s maxIterations : %d\n",message,input->maxIterations);
      sprintf(message,"%s maxIterationsPerType : %d\n",
              message,input->maxIterationsPerType);
      sprintf(message,"%s step_factor : %f\n",message,input->step_factor);
      sprintf(message,"%s sides_edge_type : %d\n",message,input->sides_edge_type);
      sprintf(message,"%s selection_edge_type : %d\n",
              message,input->selection_edge_type);
      sprintf(message,"%s countGhosts : %d\n",message,input->countGhosts);
      
      EnGPar_Log_Function(message);
    }

    //Setup the original owners arrays before balancing
    input->g->setOriginalOwners();
    unsigned int index=0;
    target_dimension = input->priorities[index];
    sd = new SDSlope;
    double tol=1.1;
    if (input->tolerances.size()>index)
      tol = input->tolerances[index];
    if (1 == PCU_Comm_Peers()) {
      printf("EnGPar ran in serial, nothing to do exiting...\n");
      return;
    }
    if (!PCU_Comm_Self() && verbosity >= 0)
      printf("Starting criteria type %d with imbalances: ",target_dimension);
    printImbalances(input->g);
    int step = 0;
    int inner_steps=0;
    double time = PCU_Time();
    double targetTime=PCU_Time();
    while (step++<input->maxIterations) {
      //runStep(tol) balances the current dimension.
      //Advance to the next dimension if the current dimesion is balanced
      //or the maximum per dimension iterations is reached.
      if (!runStep(tol) || inner_steps++ >= input->maxIterationsPerType) {
        // Set the imbalance limit for the higher priority dimension
        // while balancing lower priority dimensions to be the larger of
        // the specified imbalance (tgtMaxW) or, if it wasn't reached, the
        // current weight (maxW).
        completed_dimensions.push_back(target_dimension);
        double maxW = getMaxWeight(input->g,target_dimension);
        double tgtMaxW = getAvgWeight(input->g,target_dimension)*tol;
        maxW = ( maxW < tgtMaxW ) ? tgtMaxW : maxW;
        completed_weights.push_back(maxW);
        
        targetTime = PCU_Time()-targetTime;
        targetTime = PCU_Max_Double(targetTime);
        if (verbosity >= 0 && !PCU_Comm_Self()) {
          printf("Completed criteria type %d in %d steps and took %f seconds\n",
                 target_dimension, inner_steps, targetTime);
        }
        targetTime=PCU_Time();
        
        index++;
        if (index==input->priorities.size())
          break;
        inner_steps=0;
        target_dimension=input->priorities[index];
        delete sd;
        sd = new SDSlope;
        if (input->tolerances.size()>index)
          tol = input->tolerances[index];
        if (!PCU_Comm_Self()&&verbosity >= 0)
          printf("Starting criteria type %d with imbalances: ",target_dimension);
        printImbalances(input->g);

      }      
    }
    delete sd;
    time = PCU_Time()-time;

    if (verbosity >= 0) {
      time = PCU_Max_Double(time);
      if (!PCU_Comm_Self()) {
        if(step==input->maxIterations)
          printf("EnGPar ran to completion in %d iterations in %f seconds\n",
                 input->maxIterations, time);
        else
          printf("EnGPar converged in %d iterations in %f seconds\n",step,
                 time);
      }
    }
    if (verbosity >= 2) {
      times[1] = PCU_Max_Double(times[1]);
      distance_time = PCU_Max_Double(distance_time);
      if (!PCU_Comm_Self()) {
        printf("Migration took %f%% of the total time\n",times[1]/time*100);
        printf("Distance Computation took %f seconds\n",distance_time);
      }
    }
    if (EnGPar_Is_Log_Open())
      EnGPar_End_Function();
  }

  agi::Balancer* makeBalancer(Input* in,int v_) {
    if (EnGPar_Is_Log_Open()) {
      char message[25];
      sprintf(message,"makeBalancer\n");
      EnGPar_Log_Function(message);
      EnGPar_End_Function();
    }
    return new Balancer(in,v_,"balancer");
  }
}

