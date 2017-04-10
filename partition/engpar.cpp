#include "engpar.h"
#include <PCU.h>
namespace engpar {
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

  void evaluatePartition(agi::Ngraph* g) {

    //Vertex Imbalance
    agi::lid_t my_local = g->numLocalVtxs();
    agi::lid_t my_total = g->numTotalVtxs();
    printMaxMinAvgImb(my_local,"Local Vertices");
    printMaxMinAvgImb(my_total,"Total Vertices");
    //Edge type imbalance
  
  }
}
