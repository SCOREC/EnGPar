#include <engpar_support.h>
#include <apfGraph.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <parma.h>

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if (argc != 5 && argc != 6) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <model> <mesh> <tolerance> <render=[1:on|0:off]> [multi-edgetypes]\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }

  int isMultiEdge = 0;
  if( argc == 6 ) isMultiEdge = 1;

  apf::Mesh2* m=NULL;
  agi::Ngraph* g=NULL;
  //Load mesh
  gmi_register_mesh();
  m = apf::loadMdsMesh(argv[1],argv[2]);

  double tol = atof(argv[3]);
  const int render = atoi(argv[4]);
  if (!PCU_Comm_Self())
    printf("render vtk %s\n", render == 0 ? "off" : "on" );

  Parma_PrintPtnStats(m,"before");
  if(render)
    apf::writeVtkFiles("pre", m);

  double times[10];
  times[0] = PCU_Time();
  times[4] = PCU_Time();
  if (isMultiEdge) {
    int edges[2] = {0,2};
    //balance vtx>edge>elm
    g = agi::createAPFGraph(m,3,edges,2);
  }
  else {
    //balance vtx>elm
    g = agi::createAPFGraph(m,3,0);
  }
  times[0] = PCU_Time()-times[0];
  times[1] = PCU_Time();
  double step_factor = .1;
  engpar::Input* input = engpar::createDiffusiveInput(g,step_factor);
  input->addPriority(0,tol);
  if (isMultiEdge) {
    input->addPriority(1,tol);
  }
  input->addPriority(-1,tol);

  //Create the balancer
  int verbosity = 1;
  agi::Balancer* balancer = engpar::makeBalancer(input, verbosity);
  balancer->balance(tol);

  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;
  times[1] = PCU_Time() - times[1];

  //Ensure the graph is still valid
  agi::checkValidity(g);

  times[2] = PCU_Time();

  //Migration of original data structure
  //Only implemented for mesh
  agi::PartitionMap* map = g->getPartition();

  //map can be used to migrate the original structure
  apf::Migration* plan = new apf::Migration(m);
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
  m->migrate(plan);

  times[2] = PCU_Time() - times[2];
  times[4] = PCU_Time() - times[4];
  PCU_Max_Doubles(times,4);
  
  if (!PCU_Comm_Self()) {
    printf("Total Time: %f\n",times[4]);
    printf("Construct Time: %f\n",times[0]);
    printf("Balancing Time: %f\n",times[1]);
    printf("Repartition Time: %f\n",times[2]);
  }
  
  //Destroy graph
  agi::destroyGraph(g);

  if(render)
    apf::writeVtkFiles("post", m);
  Parma_PrintPtnStats(m,"after");
  
  //Destroy mesh
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
