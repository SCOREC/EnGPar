#include "engpar_version.h"
#include "engpar_support.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

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
FILE* old_stdout=NULL;
FILE* old_stderr=NULL;
void EnGPar_Debug_Open(std::string s) {
  char file[80];
  sprintf(file,"debug/debug%d.txt",PCU_Comm_Self());

  struct stat st = {0};
  if (stat("debug", &st) == -1) {
    mkdir("debug", 0700);
  }
  if (!PCU_Comm_Self())
    printf("Opening debug files. Output can be found in debug/debug#.txt\n");
  if (s.find("o")!=std::string::npos||s.find("e")!=std::string::npos) {
    old_stdout = stdout;
    stdout = freopen(file,"w",stdout);
    if (!stdout) {
      stdout = old_stdout;
      fprintf(stderr,"Failed to open debug files\n");
      return;
    }
  
    if (s.find("e")!=std::string::npos) {
      old_stderr = stderr;
      stderr = stdout;
      if (!stderr) {
	stderr = old_stderr;
	fprintf(stderr,"Failed to open debug files\n");
	return;
      }
    }
  }
}

void EnGPar_Debug_Close() {
  if (old_stdout)
    stdout = old_stdout;
  if (old_stderr)
    stderr = old_stderr;
}
