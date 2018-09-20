#include "engpar_metrics.h"
#include <PCU.h>
#include <engpar_support.h>

namespace engpar {
  wgt_t getWeight(agi::Ngraph* g, int dimension, bool countGhosts) {
    wgt_t w =0;
    if (dimension==-1) {
      
      //calculate the total weight of the vertices
      agi::GraphVertex* vtx;
      agi::VertexIterator* vitr = g->begin();
      while ((vtx = g->iterate(vitr)))
        w+=g->weight(vtx);
      if (countGhosts) {
      agi::GhostIterator* gitr = g->beginGhosts();
      while ((vtx = g->iterate(gitr)))
        w+=g->weight(vtx);
      }
    }
    else {
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g->begin(dimension);
      while ((edge = g->iterate(eitr))) 
        w+=g->weight(edge);
      g->destroy(eitr);
    }
    return w;
  }

  double EnGPar_Get_Imbalance(wgt_t my_val) {
    wgt_t max,total;
    double avg;
    MPI_Datatype type = MPI_DOUBLE;
    MPI_Allreduce(&my_val,&max,1,type,MPI_MAX,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&total,1,type,MPI_SUM,PCU_Get_Comm());
    avg = total*1.0/PCU_Comm_Peers();
    return max*1.0/avg;
  }

  void getImbalances(agi::Ngraph* g,char* imbalances, bool countGhosts) {
    for (agi::etype type = -1;type<g->numEdgeTypes();type++) {
      agi::wgt_t w = getWeight(g,type, countGhosts);
      double imb = EnGPar_Get_Imbalance(w);
      if (!PCU_Comm_Self())
        imbalances+=sprintf(imbalances, "%2.4f ", imb);
    }
  }
  
  void printMaxMinAvgImb(agi::lid_t my_val,std::string prefix) {
    agi::gid_t in = my_val;
    agi::gid_t max,min,total;
    double avg,imb;
    MPI_Datatype type = MPI_UNSIGNED_LONG;
    MPI_Allreduce(&in,&max,1,type,MPI_MAX,PCU_Get_Comm());
    MPI_Allreduce(&in,&min,1,type,MPI_MIN,PCU_Get_Comm());
    MPI_Allreduce(&in,&total,1,type,MPI_SUM,PCU_Get_Comm());
    avg = total*1.0/PCU_Comm_Peers();
    imb = max*1.0/avg;
    if (!PCU_Comm_Self())
      EnGPar_Status_Message("%s <max, min, avg, imb> = <%lu %lu %f %f>\n",
             prefix.c_str(),max,min,avg,imb);
  }
  void printMaxMinAvgImb(agi::wgt_t my_val,std::string prefix) {
    agi::wgt_t max,min,total;
    double avg,imb;
    MPI_Datatype type = MPI_DOUBLE;
    MPI_Allreduce(&my_val,&max,1,type,MPI_MAX,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&min,1,type,MPI_MIN,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&total,1,type,MPI_SUM,PCU_Get_Comm());
    avg = total*1.0/PCU_Comm_Peers();
    imb = max*1.0/avg;
    if (!PCU_Comm_Self())
      EnGPar_Status_Message("%s <max, min, avg, imb> = <%f %f %f %f>\n",
             prefix.c_str(),max,min,avg,imb);
  }
}
