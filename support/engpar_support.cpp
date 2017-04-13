#include "engpar_version.h"
#include "engpar_support.h"
#include <PCU.h>

bool handlingPCU=false;
void EnGPar_Initialize() {
  if (!PCU_Comm_Initialized())  {
    PCU_Comm_Init();
    handlingPCU=true;
  }
  if( !PCU_Comm_Self() )
    printf("EnGPar Git hash %s\n", engpar_version());
}

void EnGPar_Finalize() {
  if (handlingPCU &&PCU_Comm_Initialized())
    PCU_Comm_Free();
}
