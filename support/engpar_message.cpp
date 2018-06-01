#include <stdio.h>
#include <stdarg.h>
#include "engpar_support.h"
#include <sys/types.h>
#include <sys/stat.h>

FILE* old_stdout=NULL;
FILE* old_stderr=NULL;
void EnGPar_Debug_Open(std::string s) {
  if (!EnGPar_Check_Verbosity(0)) {
    return;
  }
    
  char file[80];
  sprintf(file,"debug/debug%d.txt",PCU_Comm_Self());

  struct stat st;
  if (stat("debug", &st) == -1) {
    mkdir("debug", 0700);
  }
  if (!PCU_Comm_Self())
    EnGPar_Status_Message("Opening debug files. Output can be found in debug/debug#.txt\n");
  if (s.find("o")!=std::string::npos||s.find("e")!=std::string::npos) {
    old_stdout = stdout;
    stdout = freopen(file,"w",stdout);
    if (!stdout) {
      stdout = old_stdout;
      EnGPar_Error_Message("Failed to open debug files\n");
      return;
    }
  
    if (s.find("e")!=std::string::npos) {
      old_stderr = stderr;
      stderr = stdout;
      if (!stderr) {
        stderr = old_stderr;
        EnGPar_Error_Message("Failed to open debug files\n");
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


void EnGPar_Print(FILE* f, const char* prefix, const char* format, va_list args) {
  //Dangerous if someone prints a string of size > 1024
  char buffer[1024];
  sprintf(buffer,"%s %s",prefix,format);
  vfprintf(f,buffer,args);

}

int test_flag;
void EnGPar_Start_Test() {
  test_flag=0;
}
bool EnGPar_Hit_Error() {
  return test_flag>0;
}
void EnGPar_Status_Message(const char* format, ...) {
  va_list ap;
  va_start(ap,format);
  EnGPar_Print(stdout,"ENGPAR",format,ap);
  va_end(ap);
}
void EnGPar_Status_Message(int verbosity,const char* format, ...) {
  if (!EnGPar_Check_Verbosity(verbosity))
    return;
  va_list ap;
  va_start(ap,format);
  EnGPar_Print(stdout,"ENGPAR",format,ap);
  va_end(ap);
}

void EnGPar_Warning_Message(const char* format, ...) {
  va_list ap;
  va_start(ap,format);
  EnGPar_Print(stdout,"[WARNING] ENGPAR",format,ap);
  va_end(ap);

}
void EnGPar_Error_Message(const char* format, ...) {
  test_flag=1;
  va_list ap;
  va_start(ap,format);
  EnGPar_Print(stderr,"[ERROR] ENGPAR",format,ap);
  va_end(ap);
}
