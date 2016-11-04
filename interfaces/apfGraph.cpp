#include "apfGraph.h"
#include <stdio.h>
//TODO: may want to replace with something more light weight in future
#include <vector>
#include <apfNumbering.h>
#include <map>
#include <iostream>
namespace agi {

apfGraph::apfGraph(apf::Mesh* mesh,int primary_dimension,
                   int secondary_dimension) : Ngraph() {
  checkDims(mesh->getDimension(),primary_dimension,secondary_dimension);

  m = mesh;
  setupAndNumbering(primary_dimension);

  createEdges(primary_dimension,secondary_dimension);
  
}

apfGraph::apfGraph(apf::Mesh* mesh, int primary_dimension,
                   int* secondary_dimensions, int n) {
  for (int i=0;i<n;i++) {
    checkDims(mesh->getDimension(),primary_dimension,secondary_dimensions[i]);
  }
  m=mesh;
  setupAndNumbering(primary_dimension);

  for (int i=0;i<n;i++) {
    createEdges(primary_dimension,secondary_dimensions[i]);
  }
  
}

  

void apfGraph::checkDims(int dim, int primary, int second) {
  //Error checking on primary and secondary dims
  if (primary>dim ||
      primary<0 ||
      second>dim ||
      second<0) {
    fprintf(stderr,"[ERROR] primary or secondary dimensions are invalid\n");
    throw 1;
  }
  
  if (primary==second) {
    fprintf(stderr,"[ERROR] primary and secondary dimensions cannot be equal\n");
    throw 2;
  }

}

void apfGraph::setupAndNumbering(int primary_dimension) {
  //TODO: copy over coordinates
  num_verts = m->count(primary_dimension);
  if (num_verts==0) {
    fprintf(stderr,"[ERROR] Mesh is empty, exiting\n");
    throw 3;
  }
  //Create vtx weight array
  weights = new double[num_verts];
  
  //Create a global numbering on the mesh over primary_dimension
  apf::Numbering* numbers = apf::numberOwnedDimension(m,"primary_ids",
                                                      primary_dimension);
  global_nums = apf::makeGlobal(numbers);
  apf::synchronize(global_nums);

}

  
void apfGraph::createEdges(int primary_dimension,
                           int secondary_dimension) {
  etype type = addEdgeType();
  num_edges[type] = 0;
  out_degree_list[type] = new int[num_verts+1];
  out_degree_list[type][0]=0;
  std::map<apf::MeshEntity*,int> unique_edges;
  //TODO: Find a way to conceal the Edge type
  std::vector<Edge> out_es; 
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(primary_dimension);
  while ((ent = m->iterate(itr))) {
    //Set vertex weights
    size_t vert = apf::getNumber(global_nums,ent,0);
    weights[vert]=1;
    if (!m->isOwned(ent))
      continue;
    apf::Adjacent adj;
    m->getAdjacent(ent,secondary_dimension,adj);
    int nents = adj.getSize();
    for (int i=0;i<nents;i++) {
      apf::MeshEntity* bridge = adj[i];
      apf::Adjacent neighbors;
      m->getAdjacent(bridge,primary_dimension,neighbors);
      int npri = neighbors.size();
      for (int j=0;j<npri;j++) {
        apf::MeshEntity* neighbor = neighbors[j];
        if (neighbor!=ent) {
          unique_edges[neighbor]++;
        }
      }

    }
    std::map<apf::MeshEntity*,int>::iterator itr;
    for (itr=unique_edges.begin();itr!=unique_edges.end();itr++) {
      out_es.push_back(makeEdge(vert,apf::getNumber(global_nums,itr->first,0),
                                   itr->second));
    }
    unique_edges.clear();
    if (vert<num_verts)
      out_degree_list[type][vert+1]=out_es.size();
  }
  
  num_edges[type] = out_es.size();
  out_edges[type] = new Edge[num_edges[type]];//again need to conceal the type
  std::copy(out_es.begin(),out_es.end(),out_edges[type]);  

  }


}
