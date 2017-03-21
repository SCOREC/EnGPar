
#ifndef APF_MESH
#define APF_MESH

#include "ngraph.h"
#include <apfMesh2.h>
#include <map>
namespace agi {

/** \brief Create the ngraph for a SCOREC mesh with one edge type
 * \param m the mesh
 * \param primary_dimension the mesh dimension to use for graph vertices
 * \param secondary_dimension the mesh dimension to use for graph hyperedges
 */
Ngraph* createAPFGraph(apf::Mesh* m, int primary_dimension,int secondary_dimension);
/** \brief Create the ngraph for a SCOREC mesh with multiple edge type
 * \param m the mesh
 * \param primary_dimension the mesh dimension to use for graph vertices
 * \param secondary_dimensions the mesh dimensions to use for each type of graph hyperedges
 * \param num_dimensions the number of edge types to be used (should be the size of secondary_dimensions
 */
Ngraph* createAPFGraph(apf::Mesh* m, int primary_dimension,int* secondary_dimensions,
	       int num_dimensions);
/** \brief An extension of the N-Graph for SCOREC meshes

 */
class apfGraph : public Ngraph {
 private:
  apf::Mesh* m;
  apf::GlobalNumbering* global_nums;
  apf::GlobalNumbering* edge_nums[MAX_TYPES];
  std::vector<gid_t> ghosts;
  std::vector<part_t> owns;

 public:
  // \cond INTERFACE
  apfGraph(apf::Mesh*, int primary_dimension, int secondary_dimension);
  apfGraph(apf::Mesh*, int primary_dimension, int* secondary_dimensions,int n);
  void destroyData();
  ~apfGraph();
    
  //Utility
  void migrate(std::map<GraphVertex*,int>&) {};
  // \endcond
 private:
  void checkDims(int dim,int primary,int second);
  void setupPrimary(int primary);
  etype setupSecondary(int second);
  void connectToEdges(int primary,int second, etype type);
  void connectToPins(int primary,int second, etype type);
  void constructGhostVerts();
};
 
}//agi namespace

#endif
