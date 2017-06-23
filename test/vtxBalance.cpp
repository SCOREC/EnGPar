#include <engpar_support.h>
#include <apfGraph.h>
#include <binGraph.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <engpar.h>

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if ( argc!=3 &&argc != 4&&argc!=5) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph> <step factor>\n", argv[0]);
      printf("Usage: %s <model> <mesh> <step factor> [verbosity]\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }

  //Load mesh
  agi::Ngraph* g;
  apf::Mesh2* m;
  double step_factor;
  if (argc>3) {
    gmi_register_mesh();
    m = apf::loadMdsMesh(argv[1],argv[2]);
    step_factor = atof(argv[3]);
    //Construct graph
    g = agi::createAPFGraph(m,3,0);
  }
  else {
    g = agi::createBinGraph(argv[1]);
    step_factor = atof(argv[2]);
  }
  
  engpar::evaluatePartition(g);
  
  //Create the balancer
  agi::Balancer* balancer = engpar::makeVtxBalancer(g,step_factor,argc==5);
  balancer->balance(1.01);
  engpar::evaluatePartition(g);
  //Destroy balancer
  delete balancer;
  
  //Destroy graph
  agi::destroyGraph(g);

  if (argc>3) {
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
