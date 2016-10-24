
#ifndef APF_MESH
#define APF_MESH

#include "agi.h"
#include <apfMesh2.h>
#include <map>
namespace agi {

class apfGraph : public Ngraph {
 private:
  apf::Mesh* m;
 public:
   apfGraph(apf::Mesh*, int primary_dimension, int secondary_dimension);
   ~apfGraph() {};
    
   //Vertex Operations
   double weight(GraphVertex*) const {};

   int owner(GraphVertex*) const {};

   const std::vector<double>& coordinates(GraphVertex*) const {};

   int degree(GraphVertex*,etype) const {};
    
   int* edges(GraphVertex*,etype) const {};

   //Edge Operations
   double weight(GraphEdge*) const {};

   GraphVertex* u(GraphEdge*) const {};
   GraphVertex* v(GraphEdge*) const {};

   GraphVertex* other(GraphEdge*,GraphVertex*) const {};

   //Traversal
   VertexIterator* begin() {};
   VertexIterator* iterate(VertexIterator*) {};
   void destroy(VertexIterator*) {};

   //Utility
   bool isEqual(GraphVertex*,GraphVertex*) {};
   void migrate(std::map<GraphVertex*,int>&) {};
};
 
}//agi namespace

#endif
