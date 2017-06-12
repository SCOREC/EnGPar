#include "engpar_balancer.h"
#include "engpar_queue.h"
#include <PCU.h>
#include "../engpar.h"
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
  }
  Balancer::Balancer(Input* input_, int v, const char* n) :
    agi::Balancer(input_->g,v,n), input(input_) {
    times[0]=0;
    times[1]=0;
  }
  bool Balancer::runStep(double tolerance) {
    double time[2];
    time[0] = PCU_Time();
    Sides* sides = makeSides(input);
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), sides->print("Sides").c_str());
    Weights* vtxWeights = makeVtxWeights(input, sides);
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), vtxWeights->print("Weights").c_str());
    
    Weights** edgeWeights= NULL;/* = new Weights*[input->g->numEdgeTypes()];
    for (agi::etype i=0;i<input->g->numEdgeTypes();i++) {
      edgeWeights[i] = makeEdgeWeights(sides,i);
    }
    */
    Targets* targets = makeTargets(input,sides,vtxWeights,edgeWeights);
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), targets->print("Targets").c_str());
    Queue* pq = createIterationQueue(input->g);
    Selector* selector = makeSelector(input,pq,&completed_dimensions,
                                      &completed_weights);
    agi::Migration* plan = new agi::Migration;
    wgt_t planW = 0.0;
    for (unsigned int cavSize=2;cavSize<=12;cavSize+=2) {
      planW += selector->select(targets,plan,planW,cavSize,target_dimension);
    }

    if (completed_dimensions.size()>0) {
      Midd* midd = selector->trim(targets,plan);
      selector->cancel(plan,midd);
    }
    
    time[0] = PCU_Time()-time[0];
    int numMigrate = plan->size();
    if (verbosity>=2) {
      int* counts = new int[PCU_Comm_Peers()];
      for (int i=0;i<PCU_Comm_Peers();i++)
	counts[i] = 0;
      agi::Migration::iterator itr;
      for (itr = plan->begin();itr!=plan->end();itr++)
	counts[itr->second]++;
      for (int i=0;i<PCU_Comm_Peers();i++)
	if (counts[i]>0)
	  printf("%d sending %d to %d\n",PCU_Comm_Self(),counts[i],i);
    }
    time[1] = PCU_Time();
    input->g->migrate(plan);
    time[1] = PCU_Time()-time[1];
    numMigrate = PCU_Add_Int(numMigrate);
    PCU_Max_Doubles(time,2);
    if (!PCU_Comm_Self()) {
      printf("Step took %f seconds\n",time[0]);
      printf("Migrating %d vertices took %f seconds\n",numMigrate,time[1]);
      times[0]+=time[0];
      times[1]+=time[1];
    }

    
    double imb = EnGPar_Get_Imbalance(vtxWeights->myWeight());
    //Check for completition of criteria
    delete pq;
    delete sides;
    delete vtxWeights;
    delete targets;
    delete selector;
    return imb>tolerance;
  }
  void Balancer::balance(double tol) {
    unsigned int index=0;
    target_dimension = input->priorities[index];
    if (1 == PCU_Comm_Peers()) return;

    int step = 0;
    int inner_steps=0;
    double time = PCU_Time();
    double targetTime=PCU_Time();
    while (step++<input->maxIterations) {
      if (!runStep(tol)||inner_steps++>=input->maxIterationsPerType) {


        completed_dimensions.push_back(target_dimension);
        double maxW = getMaxWeight(input->g,target_dimension);
        double tgtMaxW = getAvgWeight(input->g,target_dimension)*tol;
        maxW = ( maxW < tgtMaxW ) ? tgtMaxW : maxW;
        completed_weights.push_back(maxW);
        
	targetTime = PCU_Time()-targetTime;
	targetTime = PCU_Max_Double(targetTime);
	if (!PCU_Comm_Self()) {
	  printf("Completed dimension %d in %f seconds\n",target_dimension,
		 targetTime);
	}
	targetTime=PCU_Time();
        
        index++;
	if (index==input->priorities.size())
	  break;
        inner_steps=0;
	target_dimension=input->priorities[index];

      }
    }
    time = PCU_Time()-time;
    time = PCU_Max_Double(time);
    if (!PCU_Comm_Self()) {
      if (step==input->maxIterations)
	printf("EnGPar ran to completion in %d iterations in %f seconds\n",
	       input->maxIterations, time);
      else
	printf("EnGPar converged in %d iterations in %f seconds\n",step,
	       time);
      if (verbosity)
	printf("Migration took %f%% of the total time\n",times[1]/time*100);
    }
  }

  agi::Balancer* makeBalancer(Input* in,int v_) {
    return new Balancer(in,v_,"balancer");
  }
}

