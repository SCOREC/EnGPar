#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <cstring>
#include <apfGraph.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
bool cmpebin(char* str) {
  return strlen(str)>4&&strcmp(str+strlen(str)-5,".ebin")==0;             
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Open_Log();
  
  if ( argc!= 2 && argc!=3) {
    if ( !PCU_Comm_Self() ) {
      printf("Usage: %s <graph>\n", argv[0]);
      printf("Usage: %s <bgd_prefix>\n", argv[0]);
      printf("Usage: %s <model> <mesh>\n", argv[0]);
    }
    EnGPar_Finalize();
    assert(false);
  }
  agi::Ngraph* g=NULL;
  if (argc==3) {
    gmi_register_mesh();
    apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
    g= agi::createAPFGraph(m,3,2);
  }
  else if (cmpebin(argv[1]))
    g= agi::createBinGraph(argv[1]);
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }

  //Ensure the graph is valid
  agi::checkValidity(g);

  //Write the vtk files
  std::string filename = "cake";
  agi::writeVTK(g,filename.c_str());
  
  //Destroy graph
  agi::destroyGraph(g);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  EnGPar_Finalize();
  MPI_Finalize();
}
