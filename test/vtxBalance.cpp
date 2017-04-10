#include <apfGraph.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <engpar.h>
int main(int argc, char* argv[]) {

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  if ( argc != 3&&argc!=4) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> [verbosity]\n", argv[0]);
    MPI_Finalize();
    assert(false);
  }
  
  //Load mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  //Construct graph
  agi::Ngraph* g = agi::createAPFGraph(m,3,2);

  //Create the balancer
  agi::Balancer* balancer = engpar::makeVtxBalancer(g,1.1,argc==4);
  balancer->balance(0.1);

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

  PCU_Comm_Free();
  MPI_Finalize();
}
