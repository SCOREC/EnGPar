#ifndef WRITE_CAVITY_H
#define WRITE_CAVITY_H

#include <string>
#include <iostream>
#include <sstream>

namespace cavityWriter {
  extern std::stringstream foo;
  // get the stream handle
  void append(std::stringstream*);
  // gather the streams on a subset of processes and write to disk
  void writeToFile();
}

#endif

