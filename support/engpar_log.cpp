#include "engpar_support.h"
#include <cassert>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>


FILE* logFD=NULL;
int depth;
void EnGPar_Open_Log() {
  if (!PCU_Comm_Self())  {
    logFD = fopen("log.txt","w");
    depth=0;
  }
}
bool EnGPar_Is_Log_Open() {
  return logFD;
}

void EnGPar_Log_Function(char* message) {
  assert(depth>=0);
  std::string tab(depth*2,' ');
  fprintf(logFD,"%s%s",tab.c_str(),message);
  depth++;
}
void EnGPar_End_Function() {
  depth--;
}

void EnGPar_Close_Log() {
  if (logFD)
    fclose(logFD);
  logFD=NULL;
}
