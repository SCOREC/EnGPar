#include "ngraph.h"
#include "engpar_support.h"
#include <cstring>
#include <stdexcept>
#ifdef KOKKOS_ENABLED 
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
#endif
namespace agi {

  Ngraph* createEmptyGraph() {
    return new Ngraph;
  }
  Ngraph::Ngraph() {
    isHyperGraph=false;
    num_global_verts=0;
    num_local_verts=0;
    num_ghost_verts=0;
  
    num_types=0;
    local_weights=NULL;
    local_coords=NULL;
    for (int i=0;i<MAX_TYPES;i++) {
      edge_unmap[i] = NULL;
      edge_weights[i] = NULL;
      num_global_edges[i]=0;
      num_local_edges[i]=0;
      degree_list[i]=NULL;
      edge_list[i]=NULL;
      num_global_pins[i]=0;
      num_local_pins[i]=0;
      pin_degree_list[i]=NULL;
      pin_list[i]=NULL;
      vev_offsets[i] = NULL;
      vev_lists[i] = NULL;
      eve_offsets[i] = NULL;
      eve_lists[i] = NULL;
    }

    local_unmap = NULL;
    ghost_unmap=NULL;
    owners = NULL;
    original_owners = NULL;

    ghost_weights = NULL;
  }

  void Ngraph::constructGraph(bool isHG,
                              std::vector<gid_t>& verts,
                              std::vector<wgt_t>& wgts,
                              std::vector<gid_t>& edge_ids,
                              std::vector<lid_t>& degs,
                              std::vector<gid_t>& pins_to_verts,
                              std::unordered_map<gid_t,part_t>& owns) {
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"constructGraph: ");
      if (isHG)
        sprintf(message,"%sCreating hypergraph\n",message);
      else
        sprintf(message,"%sCreating traditional graph\n",message);
      EnGPar_Log_Function(message);
    }

    constructVerts(isHG,verts,wgts);
    std::vector<wgt_t> eWeights;
    constructEdges(edge_ids,degs,pins_to_verts,eWeights);
    constructGhosts(owns);

    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
  }
  void Ngraph::constructVerts(bool isHG,
                              std::vector<gid_t>& verts,
                              std::vector<wgt_t>& wgts) {
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"constructVerts: ");
      if (isHG)
        sprintf(message,"%sCreating hypergraph\n",message);
      else
        sprintf(message,"%sCreating traditional graph\n",message);
      EnGPar_Log_Function(message);
    }
    destroyData();
    isHyperGraph=isHG;
    num_local_verts=verts.size();
    local_unmap = new gid_t[num_local_verts];
    local_weights = new wgt_t[num_local_verts];
    local_coords = NULL;
    num_types = 0;
    for (lid_t i=0;i<(lid_t)verts.size();i++) {
      local_unmap[i] = verts[i];
      vtx_mapping[verts[i]]=i;
      //set local_weights
      if ((lid_t)wgts.size()>i)
        local_weights[i] = wgts[i];
      else
        local_weights[i] = 1;
    }

    num_ghost_verts=0;
    num_global_verts = PCU_Add_Long(num_local_verts);
    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
  }

  void Ngraph::constructVerts(bool isHG,lid_t num_verts,
                              gid_t* verts, wgt_t* wgts) {
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"constructVerts: ");
      if (isHG)
        sprintf(message,"%sCreating hypergraph\n",message);
      else
        sprintf(message,"%sCreating traditional graph\n",message);
      EnGPar_Log_Function(message);
    }
    destroyData();
    isHyperGraph=isHG;
    num_local_verts=num_verts;
    local_unmap = new gid_t[num_local_verts];
    local_weights = new wgt_t[num_local_verts];
    local_coords = NULL;
    num_types = 0;
    for (lid_t i=0;i<num_verts;i++) {
      local_unmap[i] = verts[i];
      vtx_mapping[verts[i]]=i;
      //set local_weights
      if (wgts!=NULL)
        local_weights[i] = wgts[i];
      else
        local_weights[i] = 1;
    }
    num_ghost_verts=0;
    num_global_verts = PCU_Add_Long(num_local_verts);
    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
  }

  etype Ngraph::constructEdges(std::vector<gid_t>& edge_ids,
                               std::vector<lid_t>& degs,
                               std::vector<gid_t>& pins_to_verts,
                               std::vector<wgt_t>& e_weights) {
    etype t = addEdgeType();
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"constructEdges: %d\n",t);
      EnGPar_Log_Function(message);
    }
    //Setup known array sizes
    degree_list[t] = new lid_t[num_local_verts+1];
    for (lid_t i=0;i<num_local_verts+1;i++)
      degree_list[t][i] = 0;
    num_local_edges[t] = edge_ids.size();
    num_local_pins[t] = pins_to_verts.size();
    edge_unmap[t] = new gid_t[edge_ids.size()];
    //Note pin*list is not used in traditional graph but its made here for generality.
    pin_degree_list[t] = new lid_t[degs.size()+1];
    pin_degree_list[t][0]=0;
    pin_list[t] = new lid_t[pins_to_verts.size()];
    lid_t new_ghosts=0;

    /*Loop through edges:
        determine vertex degrees
        add ghost vertices
    */
    for (lid_t i=0;i<(lid_t)edge_ids.size();i++) {
      //set edge_weight
      gid_t gid = edge_ids[i];
      edge_mapping[t][gid]=i;
      edge_unmap[t][i]=gid;
      pin_degree_list[t][i+1]=pin_degree_list[t][i]+degs[i];
      for (lid_t j=pin_degree_list[t][i];j<pin_degree_list[t][i+1];j++) {
        gid_t v = pins_to_verts[j];
        map_t::iterator vitr = vtx_mapping.find(v);
        //If the vtx is locally owned
        if (vitr!=vtx_mapping.end()&&vitr->second<num_local_verts) {
          //If in traditional mode skip the first pin (source vtx)
          if (!(j%2))
            degree_list[t][vitr->second+1]++;
          //If in hypergraph mode add to all vertices
          //   because each vtx has 1 pin to the hyperedge
          else if (isHyperGraph)
            degree_list[t][vitr->second+1]++;
          pin_list[t][j]=vitr->second;
        }
        else {//Otherwise it is a ghost
          //If the ghost is new then add it
          if (vitr==vtx_mapping.end()) {
            vtx_mapping[v]=num_local_verts+num_ghost_verts;
            pin_list[t][j] = num_local_verts+num_ghost_verts++;
            new_ghosts++;
          }
          else {
            pin_list[t][j]=vitr->second;
          }
        }
      }
    }
    //In traditional mode the number of edges is half the pins
    if (!isHyperGraph)
      num_local_edges[t] = num_local_pins[t]/2;
    //Convert degree list into a list of offsets
    for (lid_t i=1;i<num_local_verts+1;i++) {
      degree_list[t][i]+=degree_list[t][i-1];
    }
    //temp_counts copies the degree_list so we can build the edge_list incrementally
    lid_t* temp_counts = (lid_t*)malloc(num_local_verts*sizeof(lid_t));
    std::memcpy(temp_counts, degree_list[t], num_local_verts*sizeof(lid_t));
    edge_list[t] = new lid_t[degree_list[t][num_local_verts]];
    //If there are ghost vertices then we make the mapping between localid->globalid
    if (new_ghosts>0) {
      if (ghost_unmap) {
        gid_t* tmp_ghosts = new gid_t[num_ghost_verts];
        for (lid_t i=0;i<num_ghost_verts-new_ghosts;i++)
          tmp_ghosts[i]=ghost_unmap[i];
        delete [] ghost_unmap;
        ghost_unmap = tmp_ghosts;
      }
      else
        ghost_unmap = new gid_t[num_ghost_verts];
    }

    //Loop through edges again to build the edge lists
    for (lid_t i=0;i<(lid_t)edge_ids.size();i++) {
      for (lid_t j=pin_degree_list[t][i];j<pin_degree_list[t][i+1];j++) {
        lid_t u = pin_list[t][j];
        if (u>=num_local_verts){
          ghost_unmap[u-num_local_verts] = pins_to_verts[j];
        }
        if (isHyperGraph) {
          //In hypergraph mode the edges are numbered by a local id/global id
          if (u<num_local_verts)
            edge_list[t][temp_counts[u]++] = i;
        }
        else {
          //In traditional mode the edge_list stores the destination of the edge
          lid_t v = pin_list[t][++j];
          if (v>=num_local_verts){
            ghost_unmap[v-num_local_verts] = pins_to_verts[j];
          }
          if (u<num_local_verts)
            edge_list[t][temp_counts[u]++] = v;
        }
      }
    }
    free(temp_counts);

    //Set the weights of the edges
    edge_weights[t] = new wgt_t[edge_ids.size()];
    if (e_weights.size()==0)
      for (lid_t i=0;i<(lid_t)edge_ids.size();i++) 
        edge_weights[t][i] = 1;
    else
      for (lid_t i=0;i<(lid_t)e_weights.size();i++)
        edge_weights[t][i] = e_weights[i];

    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
    return t;
  }
  etype Ngraph::constructEdges(gid_t num_edges, gid_t* edge_ids,
                               lid_t* degs, gid_t* pins_to_verts,
                               wgt_t* e_weights) {
    etype t = addEdgeType();
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"constructEdges: %d\n",t);
      EnGPar_Log_Function(message);
    }
    degree_list[t] = new lid_t[num_local_verts+1];
    for (lid_t i=0;i<num_local_verts+1;i++)
      degree_list[t][i] = 0;
    num_local_edges[t] = num_edges;
    num_local_pins[t] = 0;
    for (lid_t i=0;i<num_edges;i++)
      num_local_pins[t]+=degs[i];
    edge_unmap[t] = new gid_t[num_edges];
    pin_degree_list[t] = new lid_t[num_edges+1];
    pin_degree_list[t][0]=0;
    pin_list[t] = new lid_t[num_local_pins[t]];
    lid_t new_ghosts=0;
    for (lid_t i=0;i<num_edges;i++) {
      //set edge_weight
      gid_t gid = edge_ids[i];
      edge_mapping[t][gid]=i;
      edge_unmap[t][i]=gid;
      pin_degree_list[t][i+1]=pin_degree_list[t][i]+degs[i];
      for (lid_t j=pin_degree_list[t][i];j<pin_degree_list[t][i+1];j++) {
        gid_t v = pins_to_verts[j];
        map_t::iterator vitr = vtx_mapping.find(v);
        if (vitr!=vtx_mapping.end()&&vitr->second<num_local_verts) {
          if (!(j%2))
            degree_list[t][vitr->second+1]++;
          else if (isHyperGraph)
            degree_list[t][vitr->second+1]++;
          pin_list[t][j]=vitr->second;
        }
        else {
          if (vitr==vtx_mapping.end()) {
            vtx_mapping[v]=num_local_verts+num_ghost_verts;
            pin_list[t][j] = num_local_verts+num_ghost_verts++;
            new_ghosts++;
          }
          else {
            pin_list[t][j]=vitr->second;
          }
        }
      }
    }
    if (!isHyperGraph)
      num_local_edges[t] = num_local_pins[t]/2;
    for (lid_t i=1;i<num_local_verts+1;i++) {
      degree_list[t][i]+=degree_list[t][i-1];
    }
    lid_t* temp_counts = (lid_t*)malloc(num_local_verts*sizeof(lid_t));
    std::memcpy(temp_counts, degree_list[t], num_local_verts*sizeof(lid_t));
    edge_list[t] = new lid_t[degree_list[t][num_local_verts]];
    if (new_ghosts>0) {
      if (ghost_unmap) {
        gid_t* tmp_ghosts = new gid_t[num_ghost_verts];
        for (lid_t i=0;i<num_ghost_verts-new_ghosts;i++)
          tmp_ghosts[i]=ghost_unmap[i];
        delete [] ghost_unmap;
        ghost_unmap = tmp_ghosts;
      }
      else
        ghost_unmap = new gid_t[num_ghost_verts];
    }

    for (lid_t i=0;i<num_edges;i++) {
      for (lid_t j=pin_degree_list[t][i];j<pin_degree_list[t][i+1];j++) {
        lid_t u = pin_list[t][j];
        if (u>=num_local_verts){
          ghost_unmap[u-num_local_verts] = pins_to_verts[j];
        }
        if (isHyperGraph) {
          if (u<num_local_verts)
            edge_list[t][temp_counts[u]++] = i;
        }
        else {
          lid_t v = pin_list[t][++j];
          if (v>=num_local_verts){
            ghost_unmap[v-num_local_verts] = pins_to_verts[j];
          }
          if (u<num_local_verts)
            edge_list[t][temp_counts[u]++] = v;
        }
      }
    }
    free(temp_counts);

    //Set the weights of the edges
    edge_weights[t] = new wgt_t[num_edges];
    if (!e_weights)
      for (lid_t i=0;i<num_edges;i++) 
        edge_weights[t][i] = 1;
    else
      for (lid_t i=0;i<num_edges;i++)
        edge_weights[t][i] = e_weights[i];

    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
    return t;
  }

  double weightTag(Ngraph* g, GraphVertex* v) {
    return g->weight(v);
  }

  void Ngraph::constructGhosts(std::unordered_map<gid_t,part_t>& owns) {
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"constructGhosts\n");
      EnGPar_Log_Function(message);
    }

    owners = NULL;
    if (num_ghost_verts>0)
      owners = new part_t[num_ghost_verts];
    for (lid_t v = 0;v<num_ghost_verts;v++) {
      owners[v] = owns[ghost_unmap[v]];
    }

    for (etype t = 0;t<num_types;t++)  {
      if (isHyperGraph) {
        gid_t nOwnedEdges=num_local_edges[t];
        gid_t nOwnedPins = num_local_pins[t];
        EdgeIterator* eitr = begin(t);
        GraphEdge* edge;
        while ((edge = iterate(eitr))) {
          Peers neighbors;
          getResidence(edge,neighbors);
          Peers::iterator itr;
          for (itr = neighbors.begin();itr!=neighbors.end();itr++) {
            if (*itr<PCU_Comm_Self()) {
              --nOwnedEdges;
              nOwnedPins-=degree(edge);
              break;
            }
          }
        }
        destroy(eitr);
        num_global_edges[t] = PCU_Add_Long(nOwnedEdges);
        num_global_pins[t] = PCU_Add_Long(nOwnedPins);
      }
      else {
        num_global_edges[t] = PCU_Add_Long(num_local_edges[t]);
        num_global_pins[t] = 2*num_global_edges[t];
      }
    }

    //Create the ghost weights
    if (PCU_Comm_Peers() > 1)
      ghost_weights = createDoubleGhostTag(weightTag);

    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
  }
  void Ngraph::constructGhosts(lid_t num_ghosts,gid_t* vert_ids,part_t* owns) {
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"constructGhosts\n");
      EnGPar_Log_Function(message);
    }

    owners = NULL;
    if (num_ghost_verts>0)
      owners = new part_t[num_ghost_verts];
    for (lid_t i = 0;i<num_ghosts;i++) {
      gid_t v = vert_ids[i];
      lid_t lid = vtx_mapping[v]-num_local_verts;
      if (lid<0)
        continue;
      owners[lid] = owns[i];
    }

    for (etype t = 0;t<num_types;t++)  {
      if (isHyperGraph) {
        gid_t nOwnedEdges=num_local_edges[t];
        gid_t nOwnedPins = num_local_pins[t];
        EdgeIterator* eitr = begin(t);
        GraphEdge* edge;
        while ((edge = iterate(eitr))) {
          Peers neighbors;
          getResidence(edge,neighbors);
          Peers::iterator itr;
          for (itr = neighbors.begin();itr!=neighbors.end();itr++) {
            if (*itr<PCU_Comm_Self()) {
              --nOwnedEdges;
              nOwnedPins-=degree(edge);
              break;
            }
          }
        }
        destroy(eitr);
        num_global_edges[t] = PCU_Add_Long(nOwnedEdges);
        num_global_pins[t] = PCU_Add_Long(nOwnedPins);
      }
      else {
        num_global_edges[t] = PCU_Add_Long(num_local_edges[t]);
        num_global_pins[t] = 2*num_global_edges[t];
      }
    }
    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
  }

  Ngraph::~Ngraph() {
    destroyData();
  }


  void Ngraph::removeEdges(etype t) {
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"destroyEdges %d\n",t);
      EnGPar_Log_Function(message);
    }

    num_local_edges[t] = 0;
    num_local_pins[t] = 0;
    num_global_edges[t] = 0;
    num_global_pins[t] = 0;

    delete [] edge_weights[t];
    edge_weights[t] = NULL;
    delete [] degree_list[t];
    degree_list[t] =NULL;
    delete [] edge_list[t];
    edge_list[t] = NULL;
    delete [] pin_degree_list[t];
    pin_degree_list[t] = NULL;
    delete [] pin_list[t];
    pin_list[t] = NULL;
    edge_mapping[t].clear();
    delete [] edge_unmap[t];
    edge_unmap[t] = NULL;
    if (t==num_types-1)
      num_types--;
  }
  
  void Ngraph::destroyData() {
    if (EnGPar_Is_Log_Open()) {
      char message[45];
      sprintf(message,"destroyData\n");
      EnGPar_Log_Function(message);
    }
    if (ghost_weights)
      destroyTag(ghost_weights);
    ghost_weights = NULL;
    if (local_weights) 
      delete [] local_weights;
    local_weights = NULL;
    if (local_coords)
      delete [] local_coords;
    local_coords = NULL;
    for (int i=0;i<MAX_TYPES;i++) {
      if (edge_unmap[i])
        delete [] edge_unmap[i];
      edge_unmap[i] = NULL;
      if (edge_weights[i])
        delete [] edge_weights[i];
      edge_weights[i] = NULL;
      if (degree_list[i])
        delete [] degree_list[i];
      degree_list[i] =NULL;
      if (edge_list[i])
        delete [] edge_list[i];
      edge_list[i] = NULL;
      if (pin_degree_list[i])
        delete [] pin_degree_list[i];
      pin_degree_list[i] = NULL;
      if (pin_list[i])
        delete [] pin_list[i];
      pin_list[i] = NULL;
      edge_mapping[i].clear();
      if (vev_offsets[i])
        delete [] vev_offsets[i];
      vev_offsets[i] = NULL;
      if (vev_lists[i])
        delete [] vev_lists[i];
      vev_lists[i] = NULL;
      if (eve_offsets[i])
        delete [] eve_offsets[i];
      eve_offsets[i] = NULL;
      if (eve_lists[i])
        delete [] eve_lists[i];
      eve_lists[i] = NULL;
    }
    if (local_unmap)
      delete [] local_unmap;
    local_unmap = NULL;
    if (ghost_unmap)
      delete [] ghost_unmap;
    ghost_unmap = NULL;
    if (owners)
      delete [] owners;
    owners = NULL;
    if (original_owners)
      delete [] original_owners;
    original_owners = NULL;
    vtx_mapping.clear();
    wp_map.clear();
    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
  }

  void destroyGraph(Ngraph* g) {
    if (EnGPar_Is_Log_Open()) {
      char message[25];
      sprintf(message,"destroyGraph\n");
      EnGPar_Log_Function(message);
    }

    delete g;
    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }  

  }

  void Ngraph::printStats() const {
    if (isHyperGraph)
      EnGPar_Status_Message("Hyper N-graph Status\n");
    else
      EnGPar_Status_Message("N-graph Status\n");
    EnGPar_Status_Message("Vertices: %ld\n", numGlobalVtxs());
    char edges[100];
    int n = sprintf(edges,"Edges: ");
    char* temp = edges + n;
    for (etype t = 0; t < numEdgeTypes(); t++) {
      n = sprintf(temp,"Type %d: %ld ",t,numGlobalEdges(t));
      temp+=n;
    }
    EnGPar_Status_Message("%s\n",edges);
    if (isHyperGraph) {
      int n = sprintf(edges,"Pins: ");
      char* temp = edges + n;
      for (etype t = 0; t < numEdgeTypes(); t++) {
        n = sprintf(temp,"Type %d: %ld ",t,numGlobalPins(t));
        temp+=n;
      }
      EnGPar_Status_Message("%s\n",edges);
    }
  }
  void Ngraph::makeUndirectedGraph() {
    if (isHyper())
      return;
    throw std::runtime_error("makeUndirectedGraph() not implemented yet");
  }

  void Ngraph::create_vev_adjacency(etype t, bool compress) {
    if (vev_offsets[t]) {
      delete [] vev_offsets[t];
      delete [] vev_lists[t];
    }
      
    vev_offsets[t] = new lid_t[num_local_verts+1];
    vev_offsets[t][0] = 0;
    int id = 1;
    GraphVertex* v;
    VertexIterator* vitr = begin();
    //Count the number of ajacencies and fill offset array
    while ((v = iterate(vitr))) {
      int count = 0;
      std::unordered_set<GraphVertex*> unique_adj;
      GraphIterator* gitr = adjacent(v);
      GraphVertex* other;
      while ((other = iterate(gitr))) {
        count++;
        unique_adj.insert(other);
      }
      destroy(gitr);
      unique_adj.erase(v);
      if (compress)
        vev_offsets[t][id] = vev_offsets[t][id-1] + unique_adj.size();
      else
        vev_offsets[t][id] = vev_offsets[t][id-1] + count;
      id++;
    }

    vev_lists[t] = new lid_t[vev_offsets[t][num_local_verts]];
    vitr = begin();
    id =0;
    //Populate the adjacency list by repeating the operation
    while ((v = iterate(vitr))) {
      std::unordered_set<GraphVertex*> unique_adj;
      GraphIterator* gitr = adjacent(v);
      GraphVertex* other;
      while ((other = iterate(gitr))) {
        if (compress)
          unique_adj.insert(other);
        else
          vev_lists[t][id++] = localID(other);
      }
      destroy(gitr);
      unique_adj.erase(v);
      std::unordered_set<GraphVertex*>::iterator itr;
      for (itr = unique_adj.begin(); itr != unique_adj.end(); itr++)
        vev_lists[t][id++] = localID(*itr);
    }
    
  }

  void Ngraph::create_eve_adjacency(etype t, bool compress) {
    if (eve_offsets[t]) {
      delete [] eve_offsets[t];
      delete [] eve_lists[t];
    }
      
    eve_offsets[t] = new lid_t[num_local_edges[t]+1];
    eve_offsets[t][0] = 0;
    int id = 1;
    GraphEdge* e;
    EdgeIterator* eitr = begin(t);
    //Count the number of ajacencies and fill offset array
    while ((e = iterate(eitr))) {
      int count = 0;
      std::unordered_set<GraphEdge*> unique_adj;
      GraphVertex* v;
      PinIterator* pitr = pins(e);
      while ((v = iterate(pitr))) {
        if (owner(v)!=PCU_Comm_Self())
          continue;
        GraphEdge* other;
        EdgeIterator* eitr2 = edges(v);
        while ((other = iterate(eitr2))) {
          count++;
          unique_adj.insert(other);
        }
        destroy(eitr2);
      }
      destroy(pitr);
      unique_adj.erase(e);
      if (compress)
        eve_offsets[t][id] = eve_offsets[t][id-1] + unique_adj.size();
      else
        eve_offsets[t][id] = eve_offsets[t][id-1] + count;
      id++;
    }
    destroy(eitr);

    eve_lists[t] = new lid_t[eve_offsets[t][num_local_edges[t]]];
    eitr = begin(t);
    id =0;
    //Populate the adjacency list by repeating the operation
    while ((e = iterate(eitr))) {
      std::unordered_set<GraphEdge*> unique_adj;
      GraphVertex* v;
      PinIterator* pitr = pins(e);
      while ((v = iterate(pitr))) {
        if (owner(v)!=PCU_Comm_Self())
          continue;
        GraphEdge* other;
        EdgeIterator* eitr2 = edges(v);
        while ((other = iterate(eitr2))) {
          if (compress)
            unique_adj.insert(other);
          else
            eve_lists[t][id++] = localID(other);
        }
        destroy(eitr2);
      }
      destroy(pitr);
      unique_adj.erase(e);
      std::unordered_set<GraphEdge*>::iterator itr;
      for (itr = unique_adj.begin(); itr != unique_adj.end(); itr++)
        eve_lists[t][id++] = localID(*itr);
    }
    destroy(eitr);
  }


#ifdef KOKKOS_ENABLED // {
  using engpar::LIDs;
  using engpar::hostToDevice;
  using engpar::deviceToHost;

  void Ngraph::buildSharedVtxMask(agi::lid_t numVerts, agi::lid_t numEdges,
      LIDs degree_view, LIDs edge_view,
      LIDs pin_degree_view, LIDs pin_edge_view,
      LIDs isSharedVtx) {
    LIDs isEdgeCut("isEdgeCut", numEdges);
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      int cut = 0;
      for (int i=pin_degree_view(e); i<pin_degree_view(e+1); ++i) {
        cut |= (pin_edge_view(i) >= numVerts);
      }
      isEdgeCut(e) = cut;
    });
    Kokkos::parallel_for(numVerts, KOKKOS_LAMBDA(const int v) {
      int shared = 0;
      for (int i=degree_view(v); i<degree_view(v+1); ++i) {
        shared |= isEdgeCut( edge_view(i) );
      }
      isSharedVtx(v) = shared;
    });
  }

  void Ngraph::parallel_create_eve(agi::etype t, bool boundaryOnly) {
    assert(isHyper());
    agi::PNgraph* pg = publicize();
    if (pg->eve_offsets[t]) {
      delete [] pg->eve_offsets[t];
      delete [] pg->eve_lists[t];
    }
    // Load graph info to device
    const int N = pg->num_local_edges[t];
    const int M = pg->num_local_verts;
    LIDs degree_view ("degree_view", M+1);
    LIDs edge_view ("edge_view", pg->degree_list[t][M]);
    LIDs pin_degree_view ("pin_degree_view", N+1);
    LIDs pin_edge_view ("pin_edge_view", pg->pin_degree_list[t][N]);
    hostToDevice(degree_view, pg->degree_list[t]);
    hostToDevice(edge_view, pg->edge_list[t]);
    hostToDevice(pin_degree_view, pg->pin_degree_list[t]);
    hostToDevice(pin_edge_view, pg->pin_list[t]);
    LIDs isSharedVtx("isVtxShared", M);
    if(boundaryOnly) {
      buildSharedVtxMask(M, N,
          degree_view, edge_view,
          pin_degree_view, pin_edge_view,
          isSharedVtx);
    } else {
      Kokkos::parallel_for(M, KOKKOS_LAMBDA(const int v) {
        isSharedVtx(v) = 1;
      });
    }
    // make hint for map size to avoid resize
    int numAdj = 0;
    Kokkos::parallel_reduce (M, KOKKOS_LAMBDA(const int v, int& upd) {
      for (int i=degree_view(v); i<degree_view(v+1); ++i) {
        for (int j=degree_view(v); j<degree_view(v+1); ++j) {
          if (isSharedVtx(v) && edge_view(i)!=edge_view(j))
            upd++;
        }
      } 
    }, numAdj);
    // fill map
    Kokkos::UnorderedMap<Kokkos::pair<int,int>,void> m (numAdj);
    Kokkos::parallel_for (M, KOKKOS_LAMBDA(const int v) {
      for (int i=degree_view(v); i<degree_view(v+1); ++i) {
        for (int j=degree_view(v); j<degree_view(v+1); ++j) {
          if (isSharedVtx(v) && edge_view(i)!=edge_view(j)) {
            m.insert( Kokkos::pair<int,int>(edge_view(i),edge_view(j)) );
          }
        }
      }
    });
    // create offset array
    LIDs deg ("deg", N+1);
    Kokkos::parallel_for (m.hash_capacity(), KOKKOS_LAMBDA(const int i) {
      if ( m.valid_at(i) ) {
        Kokkos::pair<int,int> p = m.key_at(i);
        Kokkos::atomic_fetch_add( &deg(p.first+1), 1);
      }
    });
    Kokkos::parallel_scan (N, KOKKOS_LAMBDA (const int i, int& upd, const bool& final) {
      const int val = deg(i+1);
      upd += val;
      if (final)
        deg(i+1) = upd;
    });
    pg->eve_offsets[t] = new int[N+1];
    deviceToHost(deg, pg->eve_offsets[t]);
    // make edge list array
    LIDs edgeList ("edgeList", pg->eve_offsets[t][N]);
    LIDs adjCount ("adjCount", N);
    Kokkos::parallel_for (m.hash_capacity(), KOKKOS_LAMBDA(const int i) {
      if ( m.valid_at(i) ) {
        Kokkos::pair<int,int> p = m.key_at(i);
        const int e = deg(p.first);
        const int idx = Kokkos::atomic_fetch_add( &adjCount(p.first), 1);
        edgeList(e+idx) = p.second;
      }
    });
    pg->eve_lists[t] = new int[pg->eve_offsets[t][N]];
    deviceToHost(edgeList, pg->eve_lists[t]);
  }
#endif // } KOKKOS_ENABLED
}
