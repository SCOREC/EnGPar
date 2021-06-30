#include <engpar_support.h>
#include "apfGraph.h"
#include <stdio.h>
#include <vector>
#include <apfNumbering.h>
#include <PCU.h>
#include <map>
#include <iostream>
namespace agi {

  Ngraph* createSerendipityGraph(apf::Mesh* m, const char* name, int order) {
    if (EnGPar_Is_Log_Open()) {
      char message[50];
      sprintf(message,"createSerendipityGraph %s with order %d\n",
              name, order);
      EnGPar_Log_Function(message);
      EnGPar_End_Function();
    }

    return new dofGraph(m, name, order);
  }

  dofGraph::dofGraph(apf::Mesh* mesh, const char* n, int ord) : apfGraph(), order(ord) {
    m = mesh;
    name = n;
    
    setupPrimary(m->getDimension());

    etype t = setupHyperedges();
    connectToEdges(t);
    connectPins(t);
    constructGhostVerts();
  }

  dofGraph::~dofGraph() {}

  etype dofGraph::setupHyperedges() {
    etype type = addEdgeType();

    //Count the number of DOF holders
    num_local_edges[type] = 0;
    num_global_edges[type] = 0;
    offset_global_edges[0] = 0;
    for (int d = 0; d < m->getDimension(); ++d) {
      if (hasDOFs(d)) {
        num_local_edges[type] += m->count(d);
        gid_t owned_count = countOwned(m,d);
        gid_t global_count = PCU_Add_Long(owned_count);
        num_global_edges[type] += PCU_Add_Long(global_count);
        offset_global_edges[d+1] = num_global_edges[type];
      }
    }

    //Create the edge array
    makeEdgeArray(type, num_local_edges[type]);

    //Create global numberings on the mesh for each DOF holder
    for (int d = 0; d < m->getDimension(); ++d) {
      if (hasDOFs(d)) {
        char buffer[200];
        sprintf(buffer, "%s_hyperedge_ids_dim%d", name, d);
        apf::Numbering* numbers = apf::numberOwnedDimension(m,buffer,d);
        edge_nums[d] = apf::makeGlobal(numbers);
        apf::synchronize(edge_nums[d]);
      }
    }

    //Create mapping between global and local ids
    //Additionally offset the global numberings on the mesh by the global number of lower dimension entities
    std::unordered_map<long, int> gids;
    lid_t lid = 0;
    for (int d = 0; d < m->getDimension(); ++d) {
      if (hasDOFs(d)) {
        apf::MeshEntity* ent;
        apf::MeshIterator* mitr = m->begin(d);
        while ((ent = m->iterate(mitr)) != NULL) {
          long num = apf::getNumber(edge_nums[d], ent, 0);
          num += offset_global_edges[d+1];
          gids[num] = d;
          apf::number(edge_nums[d], ent, 0, num);
          setEdge(lid++, num, numDOFs(d), type);
        }
        m->end(mitr);
      }
    }
    assert(lid == num_local_edges[type]);
    return type;
  }

  void dofGraph::connectToEdges(etype type) {
    //Construct the list of degree offsets
    degree_list[type] = new lid_t[num_local_verts+1];
    degree_list[type][0] = 0;

    //Ordered list of edges
    std::vector<lid_t> edges;
    apf::MeshEntity* ent;
    apf::MeshIterator* mitr = m->begin(m->getDimension());
    while ((ent = m->iterate(mitr)) != NULL) {
      if (!m->isOwned(ent))
        continue;

      //Get the global and local ids of the graph vertex
      gid_t gid = apf::getNumber(global_nums, ent, 0);
      lid_t lid = vtx_mapping[gid];

      //Get the adjacency from mesh elements to each lower dimension
      for (int d = 0; d < m->getDimension(); ++d) {
        if (hasDOFs(d)) {
          apf::Downward down;
          int ndown = m->getDownward(ent, d, down);
          for (int ind = 0; ind < ndown; ++ind) {
            apf::MeshEntity* bridge = down[ind];
            gid_t geid = apf::getNumber(edge_nums[d], bridge, 0);
            lid_t leid = edge_mapping[type][geid];
            //add to the list of edges
            edges.push_back(leid);
          }
        }
      }
      //update degree offset list
      degree_list[type][lid+1] = edges.size();
    }
    m->end(mitr);

    //Create the adjacency info from graph vertices to hyperedges
    edge_list[type] = new lid_t[edges.size()];
    std::copy(edges.begin(), edges.end(), edge_list[type]);
    num_global_pins[type] = PCU_Add_Long(edges.size());
  }

  void dofGraph::connectPins(etype type) {
    //Construct the list of edge degrees 
    lid_t* pdl = pin_degree_list[type] = new lid_t[num_local_edges[type]+1];
    lid_t nle = num_local_edges[type];
    for (int i = 0; i <= nle; ++i)
      pdl[i] = 0;

    //Communication loop to get the number of pins
    apf::MeshEntity* ent;
    apf::MeshIterator* itr;
    PCU_Comm_Begin();
    for (int d = 0; d < m->getDimension(); ++d) {
      if (hasDOFs(d)) {
        itr = m->begin(d);
        //Iterate over hyperedges for mesh dimension d
        while ((ent = m->iterate(itr)) != NULL) {
          gid_t vals[2];
          //get the global and local ids of the hyperedge
          vals[0] = apf::getNumber(edge_nums[d], ent, 0);
          lid_t lid = edge_mapping[type][vals[0]];

          //Get the number of owned adjacent graph vertices
          apf::Adjacent adj;
          m->getAdjacent(ent, m->getDimension(), adj);
          vals[1] = 0;
          for (unsigned int i = 0; i < adj.getSize(); ++i)
            vals[1] += m->isOwned(adj[i]);
          assert(pdl[lid+1] == 0);
          pdl[lid+1] = vals[1];
          //If the hyperedge is shared send the count to all remote copies
          if (m->isShared(ent)) {
            apf::Parts res;
            m->getResidence(ent, res);
            apf::Parts::iterator itr;
            for(itr = res.begin(); itr != res.end(); ++itr)
              if (*itr != PCU_Comm_Self()) {
                PCU_Comm_Pack(*itr, vals, 2*sizeof(gid_t));
              }
          }
        }
        m->end(itr);
      }
    }
    PCU_Comm_Send();

    //Receive the values on shared entities
    while (PCU_Comm_Receive()) {
      //vals[0] = gid of hyperedge
      //vals[1] = # of adjacent graph vertices
      gid_t vals[2];
      PCU_Comm_Unpack(vals, 2*sizeof(gid_t));
      //get the local id of the hyperedge
      lid_t lid = edge_mapping[type][vals[0]];
      //add to # of pins for the hyperedge
      pdl[lid+1] += vals[1];
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
    
    pin_list[type] = new lid_t[nlp];
    lid_t* pl = pin_list[type];

    //Communication loop to fill the pin list
    PCU_Comm_Begin();
    for(int d = 0; d < m->getDimension(); ++d) {
      if (hasDOFs(d)) {
        itr = m->begin(d);
        while((ent = m->iterate(itr)) != NULL) {
          //Get all the processes that have this hyperedge
          apf::Parts res;
          m->getResidence(ent, res);

          //vals[0] = gid of the hyperedge
          //vals[1] = gid of vertex for the pin
          gid_t vals[2];
          vals[0] = apf::getNumber(edge_nums[d],ent,0);
          lid_t lid = edge_mapping[type][vals[0]];

          //Get adjacency to mesh elements
          apf::Adjacent adj;
          m->getAdjacent(ent, m->getDimension(), adj);
          int nents = adj.getSize();

          //Loop over the adjacent mesh elements
          for (int i = 0; i < nents; ++i) {
            apf::MeshEntity* elm = adj[i];
            //skip if the mesh element is not owned
            if (!m->isOwned(elm))
              continue;

            //get the local and global id of the graph vertex
            vals[1] = apf::getNumber(global_nums, elm, 0);
            lid_t lvid = vtx_mapping[vals[1]];

            //add to own pin list
            pl[temp_counts[lid]++] = lvid;
            apf::Parts::iterator itr;
            //send the pin to all remote copies of the mesh entitiy
            for (itr = res.begin(); itr != res.end(); ++itr)
              if (*itr != PCU_Comm_Self()) {
                PCU_Comm_Pack(*itr, vals, 2*sizeof(gid_t));
              }
          }
        }
        m->end(itr);
      }
    }
    PCU_Comm_Send();

    //Iterate over received pins
    while (PCU_Comm_Receive()) {
      //vals[0] = gid of hyperedge
      //vals[1] = gid of vertex for pin
      gid_t vals[2];
      PCU_Comm_Unpack(vals,2*sizeof(gid_t));
      //get local id of hyperedge
      lid_t lid = edge_mapping[type][vals[0]];
      //try to find the graph vertex
      map_t::iterator itr = vtx_mapping.find(vals[1]);
      lid_t lvid;
      //If the vertex doesn't exist then it is a new ghost vertex
      if (itr == vtx_mapping.end()) {
        lvid = num_local_verts + (num_ghost_verts++);
        vtx_mapping[vals[1]] = lvid;
        ghosts.push_back(vals[1]);
        owns.push_back(PCU_Comm_Sender());
      }
      //otherwise grab the local id of the vertex
      else
        lvid = itr->second;
      //add vertex to the pin list
      pl[temp_counts[lid]++] = lvid;
    }
    
    //assert that we filled all of the expected pins
    for (lid_t i=0;i<nle;i++) {
      assert(temp_counts[i]==pdl[i+1]);
    }
    delete [] temp_counts;

  }

  bool dofGraph::hasDOFs(int dim) {
    //Conditional on finite element type (for now only serendipity)
    if (true) {
      if (dim == 0) {
        return true;
      }
      else if (dim == 1) {
        return (order >= 2);
      }
      else if (dim == 2) {
        return (order >= 4);
      }
      else if (dim == 3) {
        return false;
      }
      fprintf(stderr, "[ERROR] No support for 4 dimensional meshes!\n");
      return false;
    }
    else {
      fprintf(stderr, "[ERROR] hasDOFs not implemented for this type of finite element");
    }
    return false;
  }

  //TODO need to check if these are correct
  //TODO need to update this function for type of mesh element...
  int dofGraph::numDOFs(int dim) {
    //Conditional on finite element type (for now only serendipity)
    if (!hasDOFs(dim))
      return 0;
    if (true) {
      if (dim == 0) {
        return 1;
      }
      else if (dim == 1) {
        return order - 1;
      }
      else if (dim == 2) {
        //This is probably not correct
        return order - 3;
      }
    }
    else {
      fprintf(stderr, "[ERROR] numDOFs not implemented for this type of finite element");
    }
    return 0;
  }

}
