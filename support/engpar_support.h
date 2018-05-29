#ifndef ENGPAR_SUPPORT
#define ENGPAR_SUPPORT
#include <PCU.h>
#include <iostream>
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
void EnGPar_Status_Message(const char*, ...);
void EnGPar_Warning_Message(const char*, ...);
void EnGPar_Error_Message(const char*, ...);
#endif
