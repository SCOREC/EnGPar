#include <ngraph.h>
#include <binGraph.h>
#include <PCU.h>
#include <set>
#include <vector>
#include <Kokkos_Core.hpp>

void bfsTraditional(agi::Ngraph* g, agi::etype t=0){
  //The structures for passing over the graph
  std::set<agi::GraphVertex*> visited; //A set to hold all of the visited verticies
  std::vector<agi::GraphVertex*> curVector; //The vector of the vertecies currently being visited over
  std::vector<agi::GraphVertex*> nextVector; //The vector of the verticies adjacent to the currently visited verticies
  //Iterator for starting the graph traversal
  agi::VertexIterator* vitr = g->begin();
  agi::GraphVertex* curVtx = (g->iterate( vitr ));
  curVector.push_back(curVtx);
  while( !curVector.empty() ){
    nextVector.clear();
    std::vector<agi::GraphVertex*>::iterator curVtxIt;
    int sizeOfVector = (int) curVector.size();
    for( int i=0; i < sizeOfVector; ++i ){
      curVtx = curVector[i];
    //for( curVtxIt=curVector.begin(); curVtxIt != curVector.end(); ++curVtxIt ){
    //  curVtx = *(curVtxIt);
      std::set<agi::GraphVertex*>::iterator existVtx = visited.find(curVtx);
      if( existVtx != visited.end() ){
        continue;
      }
      agi::GraphIterator* neighbor = g->adjacent(curVtx,t);
      agi::GraphVertex* neighborVtx;
      while( neighborVtx = g->iterate(neighbor) ){
        nextVector.push_back(neighborVtx);
      }
      visited.insert(curVtx);
    }
    curVector = nextVector;
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file>",argv[0]);
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }
  //Construct Graph from binary file
  agi::Ngraph* g;

  g = agi::createBinGraph(argv[1]);

  //Run BFS algorithsm here

  bfsTraditional(g);
    

  destroyGraph(g);

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}
