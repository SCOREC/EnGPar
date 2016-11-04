
#ifndef APF_MESH
#define APF_MESH

#include "agi.h"
#include <apfMesh2.h>
#include <map>
namespace agi {

class apfGraph : public Ngraph {
 private:
  apf::Mesh* m;
  apf::GlobalNumbering* global_nums;
 public:
   apfGraph(apf::Mesh*, int primary_dimension, int secondary_dimension);
   apfGraph(apf::Mesh*, int primary_dimension, int* secondary_dimensions,int n);
   ~apfGraph() {};
    
   //Utility
   void migrate(std::map<GraphVertex*,int>&) {};

 private:
   void checkDims(int dim,int primary,int second);
   void setupAndNumbering(int primary);
   void createEdges(int primary,int second);
};
 
}//agi namespace

#endif
