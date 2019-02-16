#include "engpar_version.h"
#include "engpar_support.h"
#include <cassert>
#include <stdio.h>
#include <engpar_types.h>

bool handlingPCU=false;
void EnGPar_Initialize(int verbosity) {
  if (!PCU_Comm_Initialized())  {
    PCU_Comm_Init();
    handlingPCU=true;
  }
  EnGPar_Set_Verbosity(verbosity);
  if( !PCU_Comm_Self() && EnGPar_Check_Verbosity(0)) {
    EnGPar_Status_Message("Git hash %s\n", engpar_version());
  }

}

void EnGPar_Finalize() {
  EnGPar_Debug_Close();
  EnGPar_Close_Log();
  if (handlingPCU &&PCU_Comm_Initialized())
    PCU_Comm_Free();
}

void EnGPar_Switch_Comm(MPI_Comm comm) {
  PCU_Switch_Comm(comm);
}

int engpar_verbosity;
void EnGPar_Set_Verbosity(int v) {
  engpar_verbosity=v;
}

bool EnGPar_Check_Verbosity( int v) {
  return engpar_verbosity >= v;
}


#ifdef KOKKOS_ENABLED
namespace engpar {
  void hostToDevice(LIDs d, ENGPAR_LID_T* h) {
    LIDs::HostMirror hv = Kokkos::create_mirror_view(d);
    for (size_t i=0; i<hv.size(); ++i)
      hv(i) = h[i];
    Kokkos::deep_copy(d,hv);
  }
  void deviceToHost(LIDs d, ENGPAR_LID_T* h) {
    LIDs::HostMirror hv = Kokkos::create_mirror_view(d);
    Kokkos::deep_copy(hv,d);
    for(size_t i=0; i<hv.size(); ++i)
      h[i] = hv(i);
  }
  void degreeToOffset(CSR& c) {
    Kokkos::parallel_scan(c.off.dimension_0(),
      KOKKOS_LAMBDA (const int& e, int& upd, const bool& final) {
      const float size = c.off(e);
      if (final) c.off(e) = upd;
      upd += size;
    });
  }
  void allocateItems(CSR& c) {
    int listSize = 0;
    Kokkos::parallel_reduce(c.off.dimension_0()-1,
      KOKKOS_LAMBDA(const int e, int& size) {
        size += c.off(e);
      },
    listSize);
    std::string name = c.name + "_items";
    c.items = LIDs(name, listSize);
  }
}

#endif
