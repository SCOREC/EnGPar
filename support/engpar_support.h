#ifndef ENGPAR_SUPPORT
#define ENGPAR_SUPPORT
#include <PCU.h>
#include <iostream>
void EnGPar_Initialize();
void EnGPar_Finalize();


void EnGPar_Debug_Open(std::string s="o");
void EnGPar_Debug_Close();
#endif
