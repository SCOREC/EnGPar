#include "engpar_parmetis.h"
#include <iostream>
#include <PCU.h>
#include <unordered_map>
#include <stdexcept>
#include <ngraph.h>
#include <engpar_split_input.h>
#ifdef HAS_PARMETIS
#include <parmetis.h>
#endif
namespace engpar {
  
#ifdef HAS_PARMETIS
  //TODO: generalize to any edge type
  agi::Migration* EnGPar_ParMETIS(SplitInput* input, int target_parts, bool isLocal) {

    //Unique map of vertices for adjacency
    //  int can be used for edge weights
    typedef std::unordered_map<agi::GraphVertex*,int> VertexSet;
    
    agi::Ngraph* g = input->g;
    //Get an offset array of the vertices on each part
    idx_t* vtxdist = new idx_t[PCU_Comm_Peers()+1];
    for (int i = 0; i <= PCU_Comm_Self()+1; ++i) {
      vtxdist[i] = 0;
    }
    for (int i = PCU_Comm_Self()+1; i < PCU_Comm_Peers()+1; ++i) {
      vtxdist[i] = g->numLocalVtxs();
    }
    if (sizeof(idx_t) == sizeof(int)) 
      MPI_Allreduce(MPI_IN_PLACE, vtxdist, PCU_Comm_Peers()+1,
                    MPI_INT, MPI_SUM, PCU_Get_Comm());
    else
      MPI_Allreduce(MPI_IN_PLACE, vtxdist, PCU_Comm_Peers()+1,
                    MPI_LONG, MPI_SUM, PCU_Get_Comm());

    //Renumber the vertices based on the vtxdist array
    idx_t* gids = new idx_t[g->numLocalVtxs()];
    agi::VertexIterator* vitr = g->begin();
    agi::GraphVertex* vtx;
    int i = 0;
    PCU_Comm_Begin();
    agi::gid_t num_degs = 0;
    while((vtx = g->iterate(vitr))) {
      gids[i] = vtxdist[PCU_Comm_Self()] + i;

      //Get remotes
      agi::Peers res;
      agi::EdgeIterator* eitr = g->edges(vtx,input->edge_type);
      agi::GraphEdge* edge;
      VertexSet neighbors;      
      while ((edge = g->iterate(eitr))) {
        //Gather neighboring vertices
        agi::PinIterator* pitr = g->pins(edge);
        agi::GraphVertex* other;
        while ((other = g->iterate(pitr))) {
          if (other==vtx)
            continue;
          agi::part_t owner = g->owner(other);
          if (!isLocal || owner ==input->self)
            neighbors[other]=1;
          if (!isLocal)
            res.insert(owner);
        }
        g->destroy(pitr);
      }
      num_degs+=neighbors.size();
      g->destroy(eitr);
      res.erase(PCU_Comm_Self());
      agi::gid_t vals[2];
      vals[0] = g->globalID(vtx);
      vals[1] = gids[i];
      //Send new value of vertex to each remote copy
      agi::Peers::iterator itr;
      for (itr = res.begin();itr!=res.end();itr++) {
        PCU_Comm_Pack(*itr,vals,sizeof(agi::gid_t)*2);
      }
      ++i;
    }
    PCU_Comm_Send();
    //Set the renumbering of the ghost vertices
    std::unordered_map<agi::gid_t,agi::gid_t> ghost_gids;
    while(PCU_Comm_Receive()) {
      agi::gid_t vals[2];
      PCU_Comm_Unpack(vals,sizeof(agi::gid_t)*2);
      ghost_gids[vals[0]] = vals[1];
    }

    //Get the adjacency
    vitr=g->begin();
    idx_t* xadj = new idx_t[g->numLocalVtxs()+1];
    idx_t* adjncy = new idx_t[num_degs];
    xadj[0] = 0;
    i=0;
    agi::gid_t deg=0;
    idx_t* vwgts = new idx_t[g->numLocalVtxs()];
    idx_t* ewgts = NULL;
    while ((vtx = g->iterate(vitr))) {
      vwgts[i] = g->weight(vtx);
      agi::GraphIterator* gitr = g->adjacent(vtx,input->edge_type);
      agi::GraphVertex* other;
      VertexSet neighbors;
      while ((other = g->iterate(gitr))) {
        if (other==vtx)
          continue;
        neighbors[other]++;
      }
      for (VertexSet::iterator vsItr = neighbors.begin(); vsItr != neighbors.end(); vsItr++) {
        agi::GraphVertex* other = vsItr->first;
        if (g->owner(other)==input->self) {
          agi::lid_t lidv = g->localID(other);
          adjncy[deg++] = gids[lidv]; 
        }
        else if (!isLocal) {
          agi::gid_t gidv = g->globalID(other);
          adjncy[deg++] = ghost_gids[gidv];
        }
      }
      g->destroy(gitr);
      xadj[i + 1] = deg;
      ++i;
    }
    delete [] gids;
    printf("%ld %ld %ld",deg,num_degs,g->numLocalPins());
    assert(deg==num_degs);
    idx_t wgtflag=2;
    idx_t numflag=0;
    idx_t ncon=1;
    idx_t nparts = target_parts;
    real_t* tpwgts = new real_t[nparts];
    real_t ubvec = input->tolerance;
    for (i=0;i<nparts;i++) {
      tpwgts[i] = 1.0/nparts;
    }
    idx_t options[3];
    options[0] = 0;
    idx_t edgecut;
    idx_t* new_partition = new idx_t[g->numLocalVtxs()];
    MPI_Comm comm = PCU_Get_Comm();

    //Create Partition
    ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,vwgts,ewgts,
                         &wgtflag,&numflag,&ncon,&nparts,
                         tpwgts,&ubvec,options,&edgecut,
                         new_partition,&comm);
    delete [] vtxdist;
    delete [] xadj;
    delete [] adjncy;
    delete [] vwgts;
    delete [] tpwgts;
    //Make the migration plan
    agi::Migration* plan = new agi::Migration(g);
    i=0;
    vitr = g->begin();
    while((vtx = g->iterate(vitr))) {
      plan->insert(std::make_pair(vtx,new_partition[i++]));
    }
    delete [] new_partition;
    return plan;
  }
#else
  agi::Migration* EnGPar_ParMETIS(SplitInput*, int, bool) {
    throw std::runtime_error("ParMETIS not found\n");
  }
#endif
}
