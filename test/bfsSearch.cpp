#include <ngraph.h>
#include <binGraph.h>
#include <PCU.h>
#include <set>
#include <vector>
#include <Kokkos_Core.hpp>

void bfsTraditional(agi::Ngraph* g, agi::etype t=0){
  //The structures for passing over the graph
  std::set<agi::GraphVertex*> visited; //A set to hold all of the visited verticies
  std::vector<agi::GraphVertex*> curVtxVector; //The vector of the vertecies currently being visited over
  std::vector<agi::GraphVertex*> nextVtxVector; //The vector of the verticies adjacent to the currently visited verticies

  agi::VertexIterator* vitr = g->begin(); //Iterator for starting the graph traversal
  agi::GraphVertex* curVtx = (g->iterate(vitr)); //Gets the first vertex of the graph from the iterator
  curVtxVector.push_back(curVtx);
  while( !curVtxVector.empty() ){
    nextVtxVector.clear();
    int sizeOfVector = (int) curVtxVector.size();
    for( int i=0; i < sizeOfVector; ++i ){
      curVtx = curVtxVector[i];
      //Checks if the current vertex has been visited before
      std::set<agi::GraphVertex*>::iterator existVtx = visited.find(curVtx);
      if( existVtx != visited.end() ){
        continue;
      }

      //Adds the neighbors of the current vertex to the vector of the next group of vertices to visit 
      //and adds the current vertex to visited
      agi::GraphIterator* neighbor = g->adjacent(curVtx,t);
      agi::GraphVertex* neighborVtx;
      while( neighborVtx = g->iterate(neighbor) ){
        nextVtxVector.push_back(neighborVtx);
      }
      visited.insert(curVtx);
    }

    curVtxVector = nextVtxVector;
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
