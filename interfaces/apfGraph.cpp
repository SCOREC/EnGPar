#include <engpar_support.h>
#include "apfGraph.h"
#include <stdio.h>
#include <vector>
#include <apfNumbering.h>
#include <PCU.h>
#include <map>
#include <iostream>
namespace agi {

//TODO: make work for primary_dimension!=mesh_dimension
Ngraph* createAPFGraph(apf::Mesh* m, const char* name, int primary_dimension,
                       int secondary_dimension) {
  if (EnGPar_Is_Log_Open()) {
    char message[50];
    sprintf(message,"createAPFGraph %d with edge %d\n",
            primary_dimension,secondary_dimension);
    EnGPar_Log_Function(message);
    EnGPar_End_Function();
  }

  return new apfGraph(m, name, primary_dimension, secondary_dimension);
}
Ngraph* createAPFGraph(apf::Mesh* m, const char* name, int primary_dimension,
                       int* secondary_dimensions,int num_dimensions) {
  if (EnGPar_Is_Log_Open()) {
    char message[50];
    sprintf(message,"createAPFGraph %d with edges",primary_dimension);
    for (int i=0;i<num_dimensions;i++)
      sprintf(message,"%s%d ",message,secondary_dimensions[i]);
    sprintf(message,"%s\n",message);
    EnGPar_Log_Function(message);
    EnGPar_End_Function();
  }

  return new apfGraph(m, name, primary_dimension, secondary_dimensions, num_dimensions);
}
  
apfGraph::apfGraph(apf::Mesh* mesh, const char* n, int primary_dimension,
                   int secondary_dimension) : Ngraph() {
  global_nums=NULL;
  for (int i=0;i<MAX_TYPES;i++)
    edge_nums[i] = NULL;
  isHyperGraph=true;
  checkDims(mesh->getDimension(),primary_dimension,secondary_dimension);

  name = n;
  m = mesh;
  setupPrimary(primary_dimension);

  etype t = setupSecondary(secondary_dimension);
  connectToEdges(primary_dimension,secondary_dimension,t);
  connectToPins(primary_dimension,secondary_dimension,t);
  
  constructGhostVerts();
}

apfGraph::apfGraph(apf::Mesh* mesh, const char* name_, int primary_dimension,
                   int* secondary_dimensions, int n) : Ngraph(){
  global_nums=NULL;
  for (int i=0;i<MAX_TYPES;i++)
    edge_nums[i] = NULL;
  isHyperGraph=true;
  for (int i=0;i<n;i++) {
    checkDims(mesh->getDimension(),primary_dimension,secondary_dimensions[i]);
  }
  name = name_;
  m=mesh;
  
  setupPrimary(primary_dimension);

  for (int i=0;i<n;i++) {
    etype t = setupSecondary(secondary_dimensions[i]);
    connectToEdges(primary_dimension,secondary_dimensions[i],t);
    connectToPins(primary_dimension,secondary_dimensions[i],t);
  }
  constructGhostVerts();
}

  
apfGraph::~apfGraph() {
  if (global_nums)
    apf::destroyGlobalNumbering(global_nums);
  for (int i=0;i<MAX_TYPES;++i)
    if (edge_nums[i])
      apf::destroyGlobalNumbering(edge_nums[i]);
  
}
  
void apfGraph::checkDims(int dim, int primary, int second) {
  //Error checking on primary and secondary dims
  if (primary>dim ||
      primary<0 ||
      second>dim ||
      second<0) {
    if (!PCU_Comm_Self())
      printf("[ERROR] primary or secondary dimensions are invalid\n");
    throw 1;
  }

  if (primary==second) {
    if (!PCU_Comm_Self())
      printf("[ERROR] primary and secondary dimensions cannot be equal\n");
    throw 2;
  }
  if (primary != dim && PCU_Comm_Peers() > 1) {
    if (!PCU_Comm_Self())
      printf("[ERROR] partitioned mesh not supported for primary!=mesh->getDimension\n");
    throw 3;
  }


}

void apfGraph::setupPrimary(int primary_dimension) {
  
  num_local_verts = countOwned(m,primary_dimension);
  num_global_verts = PCU_Add_Long(num_local_verts);
  if (num_global_verts==0) {
    fprintf(stderr,"[ERROR] Mesh is empty, exiting\n");
    throw 3;
  }
  //Create vtx weight array
  local_weights = new wgt_t[num_local_verts];
  for (lid_t i=0;i<num_local_verts;i++) {
    local_weights[i]=0;
  }

  //Create vtx coordinate array
  local_coords = new coord_t[num_local_verts];

  //Create a global numbering on the mesh over primary_dimension
  char buffer[200];
  sprintf(buffer,"%s_primary_ids",name);
  apf::Numbering* numbers = apf::numberOwnedDimension(m,buffer,
                                                      primary_dimension);
  global_nums = apf::makeGlobal(numbers);
  apf::synchronize(global_nums);

  //Create mappings between global and local ids
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(primary_dimension);
  lid_t lid=0;
  local_unmap = new gid_t[num_local_verts];
  while ((ent = m->iterate(itr))) {
    if (!m->isOwned(ent))
      continue;
    gid_t gid = apf::getNumber(global_nums,ent,0);
    local_weights[lid]=1;
    apf::Vector3 vec = getLinearCentroid(m,ent);
    local_coords[lid][0] = vec.x();
    local_coords[lid][1] = vec.y();
    local_coords[lid][2] = vec.z();
    vtx_mapping[gid]=lid;
    local_unmap[lid++] = gid;
  }
  m->end(itr);
}

etype apfGraph::setupSecondary(int secondary_dimension) {
  etype type = addEdgeType();
  num_local_edges[type] = m->count(secondary_dimension);
  num_global_edges[type] = countOwned(m,secondary_dimension);
  num_global_edges[type] = PCU_Add_Long(num_global_edges[type]);

  //Create edge array
  makeEdgeArray(type,num_local_edges[type]);

  //Create a global numbering on the mesh over primary_dimension
  char buffer[200];
  sprintf(buffer, "%s_secondary_ids%d", name, type);
  apf::Numbering* numbers = apf::numberOwnedDimension(m,buffer,
                                                      secondary_dimension);
  edge_nums[type] = apf::makeGlobal(numbers);
  apf::synchronize(edge_nums[type]);

  //Create mappings between global and local ids
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(secondary_dimension);
  lid_t lid=0;
  while ((ent = m->iterate(itr))) {
    gid_t gid = apf::getNumber(edge_nums[type],ent,0);
    setEdge(lid++,gid,1.0,type);
  }
  m->end(itr);
  return type;
}

  
void apfGraph::connectToEdges(int primary_dimension,
                              int secondary_dimension,
                              etype type) {
  //Construct the list of degree offsets
  degree_list[type] = new lid_t[num_local_verts+1];
  degree_list[type][0]=0;

  //ordered list of edges
  std::vector<lid_t> edgs;
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(primary_dimension);
  //Iterate over the graph vertices
  while ((ent = m->iterate(itr))) {
    if (!m->isOwned(ent))
      continue;
    //Get the global and local ids of the graph vertex
    gid_t gid = apf::getNumber(global_nums,ent,0);
    lid_t lid = vtx_mapping[gid];
    apf::Adjacent adj;
    m->getAdjacent(ent,secondary_dimension,adj);
    int nents = adj.getSize();
    //Iterate over the adjacent hyperedges
    for (int i=0;i<nents;i++) {
      apf::MeshEntity* bridge = adj[i];
      gid_t geid = apf::getNumber(edge_nums[type],bridge,0);
      lid_t leid = edge_mapping[type][geid];
      //add to the list of edges
      edgs.push_back(leid);
    }
    //update degree offset lit
    degree_list[type][lid+1]=edgs.size();
  }
  m->end(itr);

  //create the adjacency info from vertices to edges
  edge_list[type] = new lid_t[edgs.size()];
  std::copy(edgs.begin(),edgs.end(),edge_list[type]);
  num_global_pins[type] = PCU_Add_Long(edgs.size());
}


void apfGraph::connectToPins(int primary_dimension,
                              int secondary_dimension,
                              etype type) {
  //Create the pin offset list
  lid_t* pdl = pin_degree_list[type] = new lid_t[num_local_edges[type]+1];
  lid_t nle = num_local_edges[type];
  pdl[0]=0;

  //Communication loop to get the number of pins
  PCU_Comm_Begin();
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(secondary_dimension);
  //Iterate over the hyperedges
  while ((ent = m->iterate(itr))) {
    gid_t vals[2];
    //get the global and local ids of the hyperedge
    vals[0]= apf::getNumber(edge_nums[type],ent,0);
    lid_t lid = edge_mapping[type][vals[0]];

    //Get the number of owned adjacent graph vertices
    apf::Adjacent adj;
    m->getAdjacent(ent,primary_dimension,adj);
    vals[1] =0;
    for (unsigned int i=0;i<adj.getSize();i++) {
      vals[1] += m->isOwned(adj[i]);
    }
    pdl[lid+1] = vals[1];
    //If the hyperedge is shared,
    //  then send the count to all copies of the hyperedge
    if (m->isShared(ent)) {
      apf::Parts res;
      m->getResidence(ent,res);
      apf::Parts::iterator itr;
      for (itr=res.begin();itr!=res.end();itr++)
        if (*itr!=PCU_Comm_Self()) {
          PCU_Comm_Pack(*itr,vals,2*sizeof(gid_t));
        }
    }
  }
  m->end(itr);
  PCU_Comm_Send();

  while (PCU_Comm_Receive()) {
    //vals[0] = gid of hyperedge
    //vals[1] = # of adjacent graph vertices
    gid_t vals[2];
    PCU_Comm_Unpack(vals,2*sizeof(gid_t));
    //get the local id of the hyperedge
    lid_t lid = edge_mapping[type][vals[0]];
    //add to # of pins for this hyperedge
    pdl[lid+1]+= vals[1];
  }
  //Make the pin degree list into an offset list
  for (lid_t i=2;i<=nle;i++) {
    pdl[i]+=pdl[i-1];
  }
  //copy over the pin offset list for a maleable copy
  lid_t* temp_counts = new lid_t[nle];
  std::copy(pdl,
            pdl+nle,
            temp_counts);
  lid_t nlp = num_local_pins[type] = pdl[nle];
  lid_t* pl = pin_list[type] = new lid_t[nlp];

  //Communication loop to get the pins
  PCU_Comm_Begin();
  itr= m->begin(secondary_dimension);
  //Iterate over hyperedges
  while ((ent = m->iterate(itr))) {
    //Get all the processes that have this hyperedge
    apf::Parts res;
    m->getResidence(ent,res);

    //vals[0] = gid of the hyperedge
    //vals[1] = gid of vertex for the pin
    gid_t vals[2];
    vals[0] = apf::getNumber(edge_nums[type],ent,0);
    lid_t lid = edge_mapping[type][vals[0]];
    apf::Adjacent adj;
    m->getAdjacent(ent,primary_dimension,adj);
    int nents = adj.getSize();
    //Loop over the adjacent graph vertices
    for (int i=0;i<nents;i++) {
      apf::MeshEntity* vtx = adj[i];
      //skip if the graph vertex is not owned
      if (!m->isOwned(vtx))
        continue;
      //get the local and global id of the graph vertex
      vals[1]= apf::getNumber(global_nums,vtx,0);
      lid_t lvid = vtx_mapping[vals[1]];
      //add to the pin list
      pl[temp_counts[lid]++] = lvid;
      apf::Parts::iterator itr;
      //send the pin to all copies of the hyperedge
      for (itr=res.begin();itr!=res.end();itr++)
        if (*itr!=PCU_Comm_Self()) {
          PCU_Comm_Pack(*itr,vals,2*sizeof(gid_t));
        }
    }
  }
  m->end(itr);
  PCU_Comm_Send();

  while (PCU_Comm_Receive()) {
    //vals[0] = gid of hyperedge
    //vals[1] = gid of vertex for pin
    gid_t vals[2];
    PCU_Comm_Unpack(vals,2*sizeof(gid_t));
    //get local id of hyperedge
    lid_t lid = edge_mapping[type][vals[0]];
    //try to find graph vertex
    map_t::iterator itr = vtx_mapping.find(vals[1]);
    lid_t lvid;
    //if the vertex doesn't exist, then it is a new ghost vertex
    if (itr==vtx_mapping.end()) {
      lvid = num_local_verts+(num_ghost_verts++);
      vtx_mapping[vals[1]]=lvid;
      ghosts.push_back(vals[1]);
      owns.push_back(PCU_Comm_Sender());
    }
    //otherwise we grab the local id of the vertex
    else
      lvid = itr->second;
    //add vertex to the pin list
    pl[temp_counts[lid]++]= lvid;

  }

  //assert that we filled all of the expected pins
  for (lid_t i=0;i<nle;i++) {
    assert(temp_counts[i]==pdl[i+1]);
  }
  delete [] temp_counts;
}


  double weightTagAPF(Ngraph* g, GraphVertex* v) {
    return g->weight(v);
  }

void apfGraph::constructGhostVerts() {
  //construct the ghost containers
  if (num_ghost_verts>0) {
    ghost_unmap = new gid_t[num_ghost_verts];
    std::copy(ghosts.begin(),ghosts.end(),ghost_unmap);
    owners = new part_t[num_ghost_verts];
    std::copy(owns.begin(),owns.end(),owners);
  }

  //Create the ghost weights
  if (PCU_Comm_Peers() > 1)
    ghost_weights = createDoubleGhostTag(weightTagAPF);

}
  

}
