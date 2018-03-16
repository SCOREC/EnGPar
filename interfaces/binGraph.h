#ifndef _BIN_GRAPH__
#define _BIN_GRAPH__

#include "ngraph.h"
#include <stdint.h>

namespace agi {

  /** \brief Create the N-Graph for a binary graph file
   * \param graph_file the compressed binary file name
   * \param part_file a partitioning of the vertices 
   *
   * If no part_file is provided a simple block partitioning will be applied.
   */
Ngraph* createBinGraph(char* graph_file,char* part_file =NULL);

/** \brief An extension of the N-Graph for binary graph files
 */
class binGraph : public Ngraph {
 public:
  // \cond INTERFACE
  //Construct an empty binary graph
  binGraph() : Ngraph() {}
  //Construction for a partitioned example
  // Part file contains partitioning
  // If no file is given applies a vertex block partitioning
  binGraph(char* graph_file,char* part_file=NULL);
  void destroyData();

  
  void migrate(std::map<GraphVertex*,part_t>&) {};
  void migrate(agi::EdgePartitionMap&);
  // \endcond
 private:
  //Loads edges from binary file
  etype load_edges(FILE* f, int64_t*& read_edges, int64_t& m_read);
  //Reads vertex owners from file
  int read_ranks(char* filename,int32_t* ranks);
  //Basic vertex block partitioning
  int vert_block_ranks(int32_t* ranks);
  //Exchanges the local edges to the correct processes
  int exchange_edges(int64_t m_read, int64_t* read_edges,
                     int32_t* ranks,etype t);
  //Creates the distributed csr
  int create_dist_csr(int32_t* ranks,etype t,bool createGhost = true);

  
};
 
}

#endif
