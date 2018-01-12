#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <engpar_split.h>
#include <PCU.h>
#include <sys/types.h>
#include <unistd.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfGraph.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <apfZoltan.h>
#include <apfPartition.h>
#include <pcu_util.h>

#define EDGE_TYPE 2

void switchToOriginals(int split_factor, bool& isOriginal, MPI_Comm& newComm);
MPI_Comm splitGraph(agi::Ngraph*&, apf::Mesh2*&, char* model, char* mesh, int split_factor);
void splitMesh(agi::Ngraph*&, apf::Mesh2*&, int argc, char* argv[]);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if (argc != 5) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <model> <mesh> <pumi/engpar parmetis (0/1)> <split_factor>\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }

  agi::Ngraph* g = NULL;
  apf::Mesh2* m = NULL;
  int use_engpar_parmetis = atoi(argv[3]);
  MPI_Comm newComm =0;
  double start_time = PCU_Time();
  int split_factor = atoi(argv[4]);
  if (use_engpar_parmetis) 
    newComm = splitGraph(g,m,argv[1],argv[2],split_factor);
  else
    splitMesh(g,m,argc,argv);
  double split_time = PCU_Time()-start_time;
  if (!PCU_Comm_Self())
    printf("Split from %d to %d parts in %.4f seconds\n",PCU_Comm_Peers()/split_factor,PCU_Comm_Peers(),split_time);
  //Create the input for diffusive load balancing (vtx>element)
  double step_factor = 0.1;
  engpar::Input* input_d = engpar::createDiffusiveInput(g,step_factor);
  input_d->addPriority(0,1.05);
  input_d->addPriority(-1,1.05);
  if (!PCU_Comm_Self())
    printf("\n");

  //Create and run the balancer
  agi::Balancer* balancer = engpar::makeBalancer(input_d,0);
  balancer->balance(1.1);

  if (!PCU_Comm_Self())
    printf("\nAfter Balancing\n");
  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;

  double total_time = PCU_Time()- start_time;
  total_time = PCU_Max_Double(total_time);
  if (!PCU_Comm_Self())
    printf("Total time to split and balance is: %.4f seconds\n",PCU_Comm_Peers()/split_factor,PCU_Comm_Peers(),total_time);
  
  //Ensure the graph is still valid
  agi::checkValidity(g);
  
  //Application continues:
  if (use_engpar_parmetis)
    MPI_Comm_free(&newComm);

  agi::PartitionMap* map = g->getPartition();

  //map can be used to migrate the original structure
  apf::Migration* plan = NULL;
  if (m) {
    plan = new apf::Migration(m);
    apf::GlobalNumbering* gids = m->findGlobalNumbering("primary_ids_global");
    apf::MeshIterator* mitr = m->begin(3);
    apf::MeshEntity* ent;
    while ((ent = m->iterate(mitr))) {
      agi::gid_t gid = apf::getNumber(gids,ent,0);
      assert(map->find(gid)!=map->end());
      agi::part_t target = map->find(gid)->second;
      if (target!=PCU_Comm_Self())
	plan->send(ent,target);
    }
    m->end(mitr);
    delete map;
  }
  if (use_engpar_parmetis) {
    gmi_register_mesh();
    gmi_model* model = gmi_load(argv[1]);
    m = apf::repeatMdsMesh(m,model,plan,split_factor);
    PCU_Barrier();
  }
  else
    m->migrate(plan);
  if (!PCU_Comm_Self())
    printf("\n");
  Parma_PrintPtnStats(m, "");
  
  
  agi::destroyGraph(g);
  if (m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
  
  EnGPar_Finalize();
  MPI_Finalize();
  return 0;
}


void switchToOriginals(int split_factor, bool& isOriginal, MPI_Comm& newComm) {
  int self = PCU_Comm_Self();
  int group = self%split_factor!=0;
  int groupRank = self%split_factor;
  isOriginal = self%split_factor==0;
  MPI_Comm_split(MPI_COMM_WORLD,group,groupRank,&newComm);
}


MPI_Comm splitGraph(agi::Ngraph*& g, apf::Mesh2*& m, char* model, char* mesh, int split_factor) {
  //Application code:
  bool isOriginal = false;
  MPI_Comm newComm;
  switchToOriginals(split_factor, isOriginal,newComm);

  //Calls to EnGPar:
  //Create the input (this sets up the internal communicators,
  //                    so this must be done before graph construction.)
  double tolerance = 1.1;
  agi::etype t = 0;
  engpar::Input* input_s = engpar::createSplitInput(g,newComm,MPI_COMM_WORLD, isOriginal,
                                                  split_factor,tolerance,t);
  
  if (isOriginal) {
    //Only the original parts will construct the graph
    gmi_register_mesh();
    m = apf::loadMdsMesh(model,mesh);
    g = agi::createAPFGraph(m,3,EDGE_TYPE);
    if (!PCU_Comm_Self())
      printf("\nBefore Split\n");
    engpar::evaluatePartition(g);
  }
  else {
      g = agi::createEmptyGraph();
  }

  engpar::split(input_s,engpar::GLOBAL_PARMETIS);
  delete input_s;
  if (!PCU_Comm_Self())
    printf("\nAfter Split\n");
  engpar::evaluatePartition(g);
  return newComm;
}


//This code is modified from core/test/zsplit.cc
const char* modelFile = 0;
const char* meshFile = 0;
int partitionFactor = 1;


apf::Migration* getPlan(apf::Mesh* m)
{
  apf::Splitter* splitter = apf::makeZoltanGlobalSplitter(
      m, apf::GRAPH, apf::PART_KWAY, false);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.1, partitionFactor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

void switchToOriginals()
{
  int self = PCU_Comm_Self();
  int groupRank = self / partitionFactor;
  int group = self % partitionFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void switchToAll()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}
void getConfig(int argc, char** argv)
{
  modelFile = argv[1];
  meshFile = argv[2];
  partitionFactor = atoi(argv[4]);
  PCU_ALWAYS_ASSERT(partitionFactor <= PCU_Comm_Peers());
}

void splitMesh(agi::Ngraph*& g, apf::Mesh2*& m, int argc, char* argv[]) {
  gmi_register_mesh();
  getConfig(argc,argv);
  gmi_model* model = gmi_load(modelFile);
  bool isOriginal = ((PCU_Comm_Self() % partitionFactor) == 0);
  apf::Migration* plan = 0;
  switchToOriginals();
  if (isOriginal) {
    m = apf::loadMdsMesh(model, meshFile);
    plan = getPlan(m);
  }
  switchToAll();
  m = apf::repeatMdsMesh(m, model, plan, partitionFactor);

  if (!PCU_Comm_Self())
    printf("\n");
  Parma_PrintPtnStats(m, "");
  if (!PCU_Comm_Self())
    printf("\n");  
  //Switch to graph
  g = agi::createAPFGraph(m,3,EDGE_TYPE);
  if (!PCU_Comm_Self())
    printf("\nAfter Split\n");
  engpar::evaluatePartition(g);
}
