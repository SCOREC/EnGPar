
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
  Ngraph* createAPFGraph(apf::Mesh* m, const char* name, int primary_dimension,int secondary_dimension);
/** \brief Create the ngraph for a SCOREC mesh with multiple edge type
 * \param m the mesh
 * \param primary_dimension the mesh dimension to use for graph vertices
 * \param secondary_dimensions the mesh dimensions to use for each type of graph hyperedges
 * \param num_dimensions the number of edge types to be used (should be the size of secondary_dimensions
 */
  Ngraph* createAPFGraph(apf::Mesh* m, const char* name, int primary_dimension,int* secondary_dimensions,
                         int num_dimensions);

  /** \brief Create the ngraph for a SCOREC mesh with one edge type for DOF holders given serendipity? finite elements
   * \param m the mesh
   * \param name a string to define the graph
   * \param The order of the serendipity finite elements
   */
  Ngraph* createSerendipityGraph(apf::Mesh* m, const char* name, int order);

/** \brief An extension of the N-Graph for SCOREC meshes

 */
class apfGraph : public Ngraph {
 protected:
  apf::Mesh* m;
  const char* name;
  apf::GlobalNumbering* global_nums;
  apf::GlobalNumbering* edge_nums[MAX_TYPES];
  std::vector<gid_t> ghosts;
  std::vector<part_t> owns;

  apfGraph();
 public:
  // \cond INTERFACE
  apfGraph(apf::Mesh*, const char* name, int primary_dimension, int secondary_dimension);
  apfGraph(apf::Mesh*, const char* name, int primary_dimension, int* secondary_dimensions,int n);
  ~apfGraph();
    
  //Utility
  void migrate(std::map<GraphVertex*,int>&) {};
  // \endcond
 protected:
  void checkDims(int dim,int primary,int second);
  void setupPrimary(int primary);
  etype setupSecondary(int second);
  void connectToEdges(int primary,int second, etype type);
  void connectToPins(int primary,int second, etype type);
  void constructGhostVerts();
};

class dofGraph : public apfGraph {
private:
  int order;
  gid_t offset_global_edges[4];
public:

  dofGraph(apf::Mesh*, const char* name, int ord);
  ~dofGraph();
private:
  //Finite Element Functions
  bool hasDOFs(int dim);
  int numDOFs(int dim);
  
  etype setupHyperedges();
  void connectToEdges(etype t);
  void connectPins(etype t);
};
}//agi namespace

#endif
