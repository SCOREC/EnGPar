#ifndef ENGPAR_SUPPORT
#define ENGPAR_SUPPORT
#include "engpar_types.h"
#include <PCU.h>
#include <iostream>
#ifdef KOKKOS_ENABLED
#include <Kokkos_Core.hpp>
#endif
void EnGPar_Initialize(int verbosity=0);
void EnGPar_Finalize();

void EnGPar_Switch_Comm(MPI_Comm);


void EnGPar_Debug_Open(std::string mode="o");
void EnGPar_Debug_Close();

void EnGPar_Open_Log();
bool EnGPar_Is_Log_Open();
void EnGPar_Log_Function(char*);
void EnGPar_End_Function();
void EnGPar_Close_Log();

void EnGPar_Set_Verbosity(int);
bool EnGPar_Check_Verbosity(int);

void EnGPar_Start_Test();
bool EnGPar_Hit_Error();
void EnGPar_Status_Message(const char*, ...);
void EnGPar_Status_Message(int verbosity, const char*, ...);
void EnGPar_Warning_Message(const char*, ...);
void EnGPar_Error_Message(const char*, ...);

double EnGPar_Peak_Memory();

#ifdef KOKKOS_ENABLED
namespace engpar {
typedef Kokkos::DefaultExecutionSpace exeSpace;
typedef Kokkos::View<ENGPAR_LID_T*, exeSpace::device_type> LIDs;
typedef Kokkos::View<ENGPAR_WGT_T*, exeSpace::device_type> WGTs;
/** \brief helper function to transfer a host array to a device view */
void hostToDevice(LIDs d, ENGPAR_LID_T* h);
void hostToDevice(WGTs d, ENGPAR_WGT_T* h);
/** \brief helper function to transfer a device view to a host array */
void deviceToHost(LIDs d, ENGPAR_LID_T* h);
/** \brief a list of n objects where object i has off(i+1)-off(i)
 *         related objects stored in items(off(i):off(i+1))
 */
struct csr {
  std::string name;
  int n;
  LIDs off;
  LIDs items;
  csr(std::string s, int numObjs,
      ENGPAR_LID_T* offsets, ENGPAR_LID_T* itms) {
    name=s;
    n=numObjs;
    std::string oname = name + "_offsets";
    std::string iname = name + "_items";
    off = LIDs(oname, n+1);
    items = LIDs(iname, offsets[n]);
    hostToDevice(off, offsets);
    hostToDevice(items, itms);
  }
  csr(std::string s, int numObjs) {
    name=s;
    n=numObjs;
    std::string str = name + "_offsets";
    off = LIDs(str, n+1);
  }
};
typedef struct csr CSR;
void degreeToOffset(CSR& c);
void allocateItems(CSR& c);
}
#endif

#endif
