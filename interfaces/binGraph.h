#ifndef _BIN_GRAPH__
#define _BIN_GRAPH__

#include "ngraph.h"
#include <stdint.h>

/* 
   Interface for a binary edge list file
 */
namespace agi {
  
class binGraph : public Ngraph {
 public:
  //Construct an empty binary graph
  binGraph() : Ngraph() {}
  //Construction for a serial example
  //  In parallel a vertex block partitioning occurs
  binGraph(char* graph_file);
  //Construction for a partitioned example
  // Part file contains partitioning
  binGraph(char* graph_file,char* part_file);
  ~binGraph();

  
  void migrate(std::map<GraphVertex*,part_t>&) {};
  void migrate(agi::EdgePartitionMap&);

 private:
  //Loads edges from binary file
  etype load_edges(char* filename,uint64_t*& read_edges, uint64_t& m_read);
  //Reads vertex owners from file
  int read_ranks(char* filename,int32_t* ranks);
  //Basic vertex block partitioning
  int vert_block_ranks(int32_t* ranks);
  //Exchanges the local edges to the correct processes
  int exchange_edges(uint64_t m_read, uint64_t* read_edges,
                     int32_t* ranks,etype t);
  //Creates the distributed csr
  int create_dist_csr(int32_t* ranks,etype t,bool createGhost = true);

  
};
 
}

#endif
