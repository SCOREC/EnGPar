#include "writeCavity.h"
#include "engpar_version.h"
#include <PCU.h>
#include <fstream>
#include <cassert>

namespace cavityWriter {
std::stringstream foo("");
int appendCalls = 0;

void append(std::stringstream* ss) {
  if (!appendCalls) {
    foo << "#engpar hash: " << engpar_version() << "\n";
  }
  foo << ss->str();
  appendCalls++;
}

void writeToFile() {
  int procsPerGroup = 2<<10;
  int numGroups = 1;
  if ( PCU_Comm_Peers() > procsPerGroup ) {
    numGroups = PCU_Comm_Peers()/procsPerGroup;
    assert( procsPerGroup * numGroups == PCU_Comm_Peers() );
  }
  int groupRank = PCU_Comm_Self() % procsPerGroup;
  int groupLeader = PCU_Comm_Self() - groupRank;
  std::string str = foo.str();
  int strSz = str.size();
  PCU_Comm_Begin();
  if( groupRank && strSz ) {
    PCU_COMM_PACK(groupLeader, strSz);
    PCU_Comm_Pack(groupLeader, str.c_str(), strSz);
  }
  PCU_Comm_Send();
  int strLen = 0;
  while( PCU_Comm_Listen() ) {
    int len;
    PCU_COMM_UNPACK(len);
    strLen += len;
    char *cstr = new char[len+1];
    PCU_Comm_Unpack(cstr,len);
    cstr[len]='\0';
    foo << cstr;
    delete[] cstr;
  }

  if(!groupRank) {
    std::stringstream ss;
    ss << "cavities_" << PCU_Comm_Self() << ".txt";
    std::ofstream out(ss.str());
    out << foo.str();
    out.close();
  }
}

}
