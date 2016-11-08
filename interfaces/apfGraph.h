
#ifndef APF_MESH
#define APF_MESH

#include "ngraph.h"
#include <apfMesh2.h>
#include <map>
namespace agi {

class apfGraph : public Ngraph {
 private:
  apf::Mesh* m;
  apf::GlobalNumbering* global_nums;
  apf::GlobalNumbering* edge_nums[MAX_TYPES];
  std::vector<gid_t> ghosts;
  std::vector<part_t> owns;

 public:
   apfGraph(apf::Mesh*, int primary_dimension, int secondary_dimension);
   apfGraph(apf::Mesh*, int primary_dimension, int* secondary_dimensions,int n);
   ~apfGraph();
    
   //Utility
   void migrate(std::map<GraphVertex*,int>&) {};

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
