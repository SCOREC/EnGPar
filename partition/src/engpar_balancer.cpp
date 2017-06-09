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
    agi::Balancer(g,v,n), factor(f) {
    maxStep =100;
    times[0]=0;
    times[1]=0;
  }
  bool Balancer::runStep(double tolerance) {
    double time[2];
    time[0] = PCU_Time();
    sides = makeSides();
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), sides->print("Sides").c_str());
    vtxWeights = makeVtxWeights(sides);
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), vtxWeights->print("Weights").c_str());
    edgeWeights = new Weights*[graph->numEdgeTypes()];
    for (agi::etype i=0;i<graph->numEdgeTypes();i++) {
      edgeWeights[i] = makeEdgeWeights(sides,i);
    }
    targets = makeTargets(sides,vtxWeights,edgeWeights);
    if (verbosity)
      printf("%d: %s\n",PCU_Comm_Self(), targets->print("Targets").c_str());
    Queue* pq = createIterationQueue(graph);
    selector = makeSelector(pq);
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
    graph->migrate(plan);
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
    //Check for completition of dimension
    return imb>tolerance;
  }
  void Balancer::balance(double tol) {
    target_dimension = 0;
    
    if (1 == PCU_Comm_Peers()) return;

    int step = 0;
    double time = PCU_Time();
    double targetTime=PCU_Time();
    while (step++<maxStep) {
      if (!runStep(tol)) {
	if (target_dimension==-1) {

	  completed_dimensions.push_back(target_dimension);
	  double maxW = getMaxWeight(graph,target_dimension);
	  double tgtMaxW = getAvgWeight(graph,target_dimension)*tol;
	  maxW = ( maxW < tgtMaxW ) ? tgtMaxW : maxW;
	  completed_weights.push_back(maxW);
	}
	targetTime = PCU_Time()-targetTime;
	targetTime = PCU_Max_Double(targetTime);
	if (!PCU_Comm_Self()) {
	  printf("Completed dimension %d in %f seconds\n",target_dimension,
		 targetTime);
	}
	targetTime=PCU_Time();
	if (target_dimension<0)
	  break;
	target_dimension--;

      }
    }
    time = PCU_Time()-time;
    time = PCU_Max_Double(time);
    if (!PCU_Comm_Self()) {
      if (step==maxStep)
	printf("EnGPar ran to completion in %d iterations in %f seconds\n",
	       maxStep, time);
      else
	printf("EnGPar converged in %d iterations in %f seconds\n",step,
	       time);
      if (verbosity)
	printf("Migration took %f%% of the total time\n",times[1]/time*100);
    }
  }
}
