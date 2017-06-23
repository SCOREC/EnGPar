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

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if ( argc!= 2 && argc != 3&& argc!=4) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph>\n", argv[0]);
      printf("Usage: %s <model> <mesh> [verbosity]\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }
  apf::Mesh2* m;
  agi::Ngraph* g;
  if (argc>2) {
    //Load mesh
    gmi_register_mesh();
    m = apf::loadMdsMesh(argv[1],argv[2]);

    //visualize the mesh before balancing
    apf::writeVtkFiles("pre", m);

    //Construct graph
    g = agi::createAPFGraph(m,3,0);
  }
  else
    g= agi::createBinGraph(argv[1]);
     
  engpar::evaluatePartition(g);

  engpar::Input* input = new engpar::Input(g);
  input->priorities.push_back(0);
  input->tolerances.push_back(1.1);
  input->priorities.push_back(-1);
  input->tolerances.push_back(1.1);
  input->step_factor=.2;
  
  //Create the balancer
  agi::Balancer* balancer = engpar::makeBalancer(input,argc==4);
  balancer->balance(1.1);

  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;

  //Migration of original data structure
  //Only implemented for mesh
  if (argc>2) {
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
  }
  
  //Destroy graph
  agi::destroyGraph(g);

  if (argc>2) {
    //visualize the mesh after balancing
    apf::writeVtkFiles("post", m);

    //Destroy mesh
    m->destroyNative();
    apf::destroyMesh(m);

  }
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
