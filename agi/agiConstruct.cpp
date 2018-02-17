#include "ngraph.h"
#include <engpar_support.h>
#include <cstring>
namespace agi {

  Ngraph* createEmptyGraph() {
    return new Ngraph;
  }
  Ngraph::Ngraph() {
    isHyperGraph=false;
    isSellCSigma=false;
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
    }

    local_unmap = NULL;
    ghost_unmap=NULL;
    owners = NULL;
    original_owners = NULL;
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
    degree_list[t] = new lid_t[num_local_verts+1];
    for (lid_t i=0;i<num_local_verts+1;i++)
      degree_list[t][i] = 0;
    num_local_edges[t] = edge_ids.size();
    num_local_pins[t] = pins_to_verts.size();
    edge_unmap[t] = new gid_t[edge_ids.size()];
    pin_degree_list[t] = new lid_t[degs.size()+1];
    pin_degree_list[t][0]=0;
    pin_list[t] = new lid_t[pins_to_verts.size()];
    lid_t new_ghosts=0;
    for (lid_t i=0;i<(lid_t)edge_ids.size();i++) {
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
    uint64_t* temp_counts = (uint64_t*)malloc(num_local_verts*sizeof(uint64_t));
    std::memcpy(temp_counts, degree_list[t], num_local_verts*sizeof(uint64_t));
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

    for (lid_t i=0;i<(lid_t)edge_ids.size();i++) {
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
    uint64_t* temp_counts = (uint64_t*)malloc(num_local_verts*sizeof(uint64_t));
    std::memcpy(temp_counts, degree_list[t], num_local_verts*sizeof(uint64_t));
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
}
