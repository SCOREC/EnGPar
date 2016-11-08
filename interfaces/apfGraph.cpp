#include "apfGraph.h"
#include <stdio.h>
//TODO: may want to replace with something more light weight in future
#include <vector>
#include <apfNumbering.h>
#include <PCU.h>
#include <map>
#include <iostream>
namespace agi {

//TODO: make work for primary_dimension!=mesh_dimension

apfGraph::apfGraph(apf::Mesh* mesh,int primary_dimension,
                   int secondary_dimension) : Ngraph() {
  primary_dimension = mesh->getDimension();
  checkDims(mesh->getDimension(),primary_dimension,secondary_dimension);

  m = mesh;
  setupPrimary(primary_dimension);

  etype t = setupSecondary(secondary_dimension);
  connectToEdges(primary_dimension,secondary_dimension,t);
  connectToPins(primary_dimension,secondary_dimension,t);

  //createGhostLayer(primary_dimension,secondary_dimension);
  
}

apfGraph::apfGraph(apf::Mesh* mesh, int primary_dimension,
                   int* secondary_dimensions, int n) {
  primary_dimension = mesh->getDimension();
  for (int i=0;i<n;i++) {
    checkDims(mesh->getDimension(),primary_dimension,secondary_dimensions[i]);
  }
  m=mesh;
  
  setupPrimary(primary_dimension);

  for (int i=0;i<n;i++) {
    etype t = setupSecondary(secondary_dimensions[i]);
    connectToEdges(primary_dimension,secondary_dimensions[i],t);
    connectToPins(primary_dimension,secondary_dimensions[i],t);
  }
  
}

apfGraph::~apfGraph() {
  //TODO: destory global numbering over primary dim and secondary dims
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

void apfGraph::setupPrimary(int primary_dimension) {
  //TODO: copy over coordinates
  
  num_local_verts = countOwned(m,primary_dimension);
  num_global_verts = num_local_verts;
  //PCU_Max_Long(num_global_verts);
  if (num_global_verts==0) {
    fprintf(stderr,"[ERROR] Mesh is empty, exiting\n");
    throw 3;
  }
  //Create vtx weight array
  local_weights = new double[num_local_verts];
  for (int i=0;i<num_local_verts;i++) {
    local_weights[i]=0;
  }

  //Create a global numbering on the mesh over primary_dimension
  apf::Numbering* numbers = apf::numberOwnedDimension(m,"primary_ids",
                                                      primary_dimension);
  global_nums = apf::makeGlobal(numbers);
  apf::synchronize(global_nums);

  //Create mappings between global and local ids
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(primary_dimension);
  lid_t lid=0;
  local_unmap = new gid_t[num_local_verts];
  while ((ent = m->iterate(itr))) {
    //Set vertex weights
    if (!m->isOwned(ent))
      continue;
    gid_t gid = apf::getNumber(global_nums,ent,0);
    local_weights[lid]=1;
    vtx_mapping[gid]=lid;
    local_unmap[lid++] = gid;
  }
}

etype apfGraph::setupSecondary(int secondary_dimension) {
  etype type = addEdgeType();
  num_local_edges[type] = m->count(secondary_dimension);
  num_global_edges[type] = countOwned(m,secondary_dimension);
  //PCU_Max_Long(num_global_edges[type]);

  //Create edge array
  makeEdgeArray(type,num_local_edges[type]);

  //Create a global numbering on the mesh over primary_dimension
  char name[20];
  sprintf(name,"secondary_ids%d",type);
  apf::Numbering* numbers = apf::numberOwnedDimension(m,name,
                                                      secondary_dimension);
  edge_nums[type] = apf::makeGlobal(numbers);
  apf::synchronize(edge_nums[type]);

  //Create mappings between global and local ids
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(secondary_dimension);
  lid_t lid=0;
  while ((ent = m->iterate(itr))) {
    //Set edge weights
    gid_t gid = apf::getNumber(edge_nums[type],ent,0);
    setEdge(lid++,gid,1.0,type);
  }
  return type;
}

  
void apfGraph::connectToEdges(int primary_dimension,
                              int secondary_dimension,
                              etype type) {
  degree_list[type] = new lid_t[num_local_verts+1];
  degree_list[type][0]=0;
  std::vector<lid_t> edgs;
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(primary_dimension);
  while ((ent = m->iterate(itr))) {
    if (!m->isOwned(ent))
      continue;
    gid_t gid = apf::getNumber(global_nums,ent,0);
    lid_t lid = vtx_mapping[gid];
    apf::Adjacent adj;
    m->getAdjacent(ent,secondary_dimension,adj);
    int nents = adj.getSize();
    for (int i=0;i<nents;i++) {
      apf::MeshEntity* bridge = adj[i];
      gid_t geid = apf::getNumber(edge_nums[type],bridge,0);
      lid_t leid = edge_mapping[type][geid];
      edgs.push_back(leid);
    }
    if (lid<num_local_verts)
      degree_list[type][lid+1]=edgs.size();
  }
  

  edge_list[type] = new lid_t[edgs.size()];
  std::copy(edgs.begin(),edgs.end(),edge_list[type]);

  
}


void apfGraph::connectToPins(int primary_dimension,
                              int secondary_dimension,
                              etype type) {

  pin_degree_list[type] = new lid_t[num_local_edges[type]+1];
  pin_degree_list[type][0]=0;
  std::vector<lid_t> edgs;
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(secondary_dimension);
  while ((ent = m->iterate(itr))) {
    if (!m->isOwned(ent))
      continue;
    gid_t gid = apf::getNumber(edge_nums[type],ent,0);
    lid_t lid = edge_mapping[type][gid];
    apf::Adjacent adj;
    m->getAdjacent(ent,primary_dimension,adj);
    int nents = adj.getSize();
    for (int i=0;i<nents;i++) {
      apf::MeshEntity* vtx = adj[i];
      gid_t gvid = apf::getNumber(global_nums,vtx,0);
      lid_t lvid = vtx_mapping[gvid];
      edgs.push_back(lvid);
    }
    if (lid<num_local_edges[type])
      pin_degree_list[type][lid+1]=edgs.size();
  }
  

  pin_list[type] = new lid_t[edgs.size()];
  std::copy(edgs.begin(),edgs.end(),pin_list[type]);
}

// void apfGraph::createGhostLayer(int primary, int secondary) {
  
//   PCU_Comm_Begin();
//   apf::MeshEntity* ent;
//   apf::MeshIterator* itr= m->begin(secondary_dimension);
//   int dim = m->getDimension();
//   while ((ent = m->iterate(itr))) {
//     if (m->isShared(ent)) {
//       apf::Copies remotes;
//       m->getRemotes(down[i],remotes);
    
//       gid_t vals[2];
//       if (dim==primary) {
//         apf::Adjacent adj;
//         m->getAdjacent(ent,primary,adj);
//         for (int i=0;i<adj.size();i++) {
//           vals[0] = apf::getNumber(global_nums,adj[i],0);
//           apf::Copies::iterator itr;
//           for (itr = remotes.begin();itr!=remotes.end();itr++) {
//             vals[1] = (uintptr_t)itr->second;
//             PCU_Comm_Pack(itr->first,vals,2);
//           }
          
//         }
//       } 
//       else {
//         //TODO: ghosting when primary!=mesh dim
//       }
//     }
//     else if (secondary==dim) {
//       //TODO: ghosting when secondary=mesh_dim
//     }
//   }
//   PCU_Comm_Send();

//   //TODO: unordered?
//   typedef std::map<apf::MeshEntity*,int> EdgeW;
//   EdgeW* unique_edges = new EdgeW[num_local_verts];
//   std::vector<part_t> ghost_owners;
//   while (PCU_Comm_Receive()) {
//     gid_t vals[2];
//     PCU_Comm_Unpack(vals,2);
//     std::map<gid_t,lid_t>::iterator itr = global_local_mapping.find(vals[0]);
//     if (itr=global_local_mapping.end()) {
//       itr = global_local_mapping.insert(std::make_pair(vals[0],
//                                                        num_local_verts+
//                                                        (num_ghost_verts++)));
//       ghost_owners.push_back(PCU_Comm_Sender());
//     }
//     lid_t ghost_lid = itr->second;
//     MeshEntity* second = (MeshEntity*)vals[1];
//     gid_t gid = apf::getNumber(global_nums,second,0);
//     unique_edges[global_local_mapping[gid]][ghost_lid]++;
//   }
//   delete [] unique_edges;

//   //TODO: get ghost weights?
  
// }
  

}
