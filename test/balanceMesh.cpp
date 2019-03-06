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
#include <engpar_diffusive_input.h>

void setWeights(agi::Ngraph* g, agi::etype t) {
  if(!PCU_Comm_Self()) {
    agi::lid_t n = g->numLocalVtxs();
    agi::wgt_t* w = new agi::wgt_t[n];
    for(int i=0; i<n; i++)
      w[i] = 2;
    g->setWeights(w);
    delete [] w;

    agi::PNgraph* pg = g->publicize();
    agi::lid_t m = g->numLocalEdges(t);
    for(int i=0; i<m; i++)
      pg->edge_weights[t][i] = 2;
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  Kokkos::initialize(argc,argv);
  EnGPar_Open_Log();
  
  if (argc != 7 && argc != 8) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <model> <mesh> <tolerance> <render=[1:on|0:off]> "
          "<kkSelect=[1:on|0:off]> <skewWeights=[1:on|0:off]> "
          "[multi-edgetypes]\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }

  int kkselect = (atoi(argv[5]) > 0);
  int skewWeights = (atoi(argv[6]) > 0);
  int isMultiEdge = 0;
  if( argc == 8 ) isMultiEdge = 1;

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
  std::string name;
  if (isMultiEdge) {
    name = "edge_02";
    int edges[2] = {0,2};
    //balance vtx>edge>elm
    g = agi::createAPFGraph(m,name.c_str(),m->getDimension(),edges,2);
  }
  else {
    name = "edge_0";
    //balance vtx>elm
    g = agi::createAPFGraph(m,name.c_str(),m->getDimension(),0);
  }
  if( skewWeights )
    setWeights(g,0);
  times[0] = PCU_Time()-times[0];
  times[1] = PCU_Time();
  double step_factor = .1;
  engpar::DiffusiveInput* input = engpar::createDiffusiveInput(g,step_factor);
  if (isMultiEdge) {
    input->addPriority(1,tol);
  }
  input->addPriority(-1,tol);
  input->maxIterationsPerType=0;
  input->kkSelect=kkselect;

  engpar::evaluatePartition(g);
  EnGPar_Debug_Open();
  //Create the balancer
  int verbosity = 1;
  engpar::balance(input,verbosity);

  engpar::evaluatePartition(g);
  times[1] = PCU_Time() - times[1];

  //Ensure the graph is still valid
  agi::checkValidity(g);

  times[2] = PCU_Time();

  //Migration of original data structure
  //Only implemented for mesh
  agi::PartitionMap* map = g->getPartition();

  //map can be used to migrate the original structure
  apf::Migration* plan = new apf::Migration(m);
  std::string numberingName = name + "_primary_ids_global";
  apf::GlobalNumbering* gids = m->findGlobalNumbering(numberingName.c_str());
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
    printf("Total Time (s): %f\n",times[4]);
    printf("Construct N-Graph from Mesh Time (s): %f\n",times[0]);
    printf("EnGPar Balancing Time (s): %f\n",times[1]);
    printf("Mesh Repartition Time (s): %f\n",times[2]);
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

  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();
}
