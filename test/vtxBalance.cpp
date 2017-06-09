#include <engpar_support.h>
#include <apfGraph.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <engpar.h>

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  if ( argc != 4&&argc!=5) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <step factor> [verbosity]\n", argv[0]);
    EnGPar_Finalize();
    assert(false);
  }

  //Load mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  //Construct graph
  agi::Ngraph* g = agi::createAPFGraph(m,3,0);

  engpar::evaluatePartition(g);
  
  //Create the balancer
  agi::Balancer* balancer = engpar::makeVtxBalancer(g,atof(argv[3]),argc==5);
  balancer->balance(1.1);

  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;
  
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
