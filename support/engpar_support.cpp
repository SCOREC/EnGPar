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
  void hostToDevice(LIDs d, AGI_LID_T* h) {
    LIDs::HostMirror hv = Kokkos::create_mirror_view(d);
    for (size_t i=0; i<hv.size(); ++i)
      hv(i) = h[i];
    Kokkos::deep_copy(d,hv);
  }
  void deviceToHost(LIDs d, AGI_LID_T* h) {
    LIDs::HostMirror hv = Kokkos::create_mirror_view(d);
    Kokkos::deep_copy(hv,d);
    for(size_t i=0; i<hv.size(); ++i)
      h[i] = hv(i);
  }
}

#endif
