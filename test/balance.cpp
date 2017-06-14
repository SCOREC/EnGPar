#include <engpar_support.h>
#include <apfGraph.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <engpar.h>
#include <engpar_input.h>
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if ( argc != 3&& argc!=4) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> [verbosity]\n", argv[0]);
    EnGPar_Finalize();
    assert(false);
  }

  //Load mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  //Construct graph
  agi::Ngraph* g = agi::createAPFGraph(m,3,0);

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

  agi::PartitionMap* map = g->getPartition();
  //map can be used to migrate the original structure
  delete map;
  
  //Destroy graph
  agi::destroyGraph(g);
  //Destroy mesh
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
