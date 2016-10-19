
#ifndef APF_MESH
#define APF_MESH

#include "agi.h"
#include <apfMesh.h>
#include <map>
namespace agi {

class APFGraph : public Ngraph {
 private:
  apf::Mesh* m;
 public:
  APFGraph(apf::Mesh*, int primary_dimension, int secondary_dimension);
    ~APFGraph();

    //Part Information
    
    int numVtx() const;
    int numEdges() const;

    
    //Vertex Operations
    virtual double weight(GraphVertex*) const;

    virtual int owner(GraphVertex*) const;

    virtual const std::vector<double>& coordinates(GraphVertex*) const;

    virtual int degree(GraphVertex*,etype) const;
    
    virtual int* edges(GraphVertex*,etype) const;

    //Edge Operations
    virtual double weight(GraphEdge*) const;

    virtual GraphVertex* u(GraphEdge*) const;
    virtual GraphVertex* v(GraphEdge*) const;

    virtual GraphVertex* other(GraphEdge*,GraphVertex*) const;

    //Traversal
    virtual VertexIterator* begin();
    virtual VertexIterator* iterate(VertexIterator*);
    virtual void destroy(VertexIterator*);

    //Utility
    virtual bool isEqual(GraphVertex*,GraphVertex*);
    virtual void migrate(std::map<GraphVertex*,int>&);
};
 
}//agi namespace

#endif
