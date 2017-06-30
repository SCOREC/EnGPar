#include "engpar.h"
#include <PCU.h>
namespace engpar {
  wgt_t getWeight(agi::Ngraph* g, int dimension) {
    wgt_t w =0;
    if (dimension==-1) {
      //calculate the total weight of the vertices
      agi::GraphVertex* vtx;
      agi::VertexIterator* vitr = g->begin();
      while ((vtx = g->iterate(vitr))) 
        w+=g->weight(vtx);
    }
    else {
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g->begin(dimension);
      while ((edge = g->iterate(eitr))) 
        w+=g->weight(edge);
      
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
  
  void printMaxMinAvgImb(agi::gid_t my_val,std::string prefix) {
    agi::gid_t max,min,total;
    double avg,imb;
    MPI_Datatype type = MPI_UNSIGNED_LONG;
    MPI_Allreduce(&my_val,&max,1,type,MPI_MAX,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&min,1,type,MPI_MIN,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&total,1,type,MPI_SUM,PCU_Get_Comm());
    avg = total*1.0/PCU_Comm_Peers();
    imb = max*1.0/avg;
    if (!PCU_Comm_Self())
      printf("%s <max, min, avg, imb> = <%lu %lu %f %f>\n",
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
      printf("%s <max, min, avg, imb> = <%f %f %f %f>\n",
             prefix.c_str(),max,min,avg,imb);
  }

  void evaluatePartition(agi::Ngraph* g) {

    //Vertex Imbalance
    agi::wgt_t my_local = getWeight(g,-1);
    agi::lid_t my_total = g->numTotalVtxs();
    printMaxMinAvgImb(my_local,"Local Vertices");
    printMaxMinAvgImb(my_total,"Total Vertices");
    //Edge type imbalance
    for (agi::etype t = 0;t<g->numEdgeTypes();t++) {
      agi::wgt_t my_edges = getWeight(g,t);
      char edge_name[15];
      sprintf(edge_name,"Local Edges %d",t);
      printMaxMinAvgImb(my_edges,edge_name);
    }
  
  }
}
