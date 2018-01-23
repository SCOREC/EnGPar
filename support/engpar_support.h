#ifndef ENGPAR_SUPPORT
#define ENGPAR_SUPPORT
#include <PCU.h>
#include <iostream>
void EnGPar_Initialize();
void EnGPar_Finalize();

void EnGPar_Switch_Comm(MPI_Comm);


void EnGPar_Debug_Open(std::string mode="o");
void EnGPar_Debug_Close();

void EnGPar_Open_Log();
bool EnGPar_Is_Log_Open();
void EnGPar_Log_Function(char*);
void EnGPar_End_Function();
void EnGPar_Close_Log();

#endif
