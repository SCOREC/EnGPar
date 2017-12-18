#include "agiMigrationTimers.h"
#include <PCU.h>
#include <cassert>

namespace agi {
  MigrationTimers::MigrationTimers() {
    int initSize = 16;
    times.reserve(initSize);
    counts.reserve(initSize);
    names.reserve(initSize);
  }

  MigrationTimers::~MigrationTimers() {
  }

  void MigrationTimers::addTimer(std::string name) {
    int idx = nameToIdx.size()+1;
    nameToIdx[name] = idx;
    names.push_back(name);
    counts.push_back(0);
    times.push_back(0);
  }

  void MigrationTimers::update(std::string name, double val) {
    int idx = getTimerIdx(name);
    times[idx] += val;
    counts[idx]++;
  }

  double MigrationTimers::processMax(std::string name) {
    int idx = getTimerIdx(name);
    return PCU_Max_Double(times[idx]);
  }

  double MigrationTimers::perCallProcessMax(std::string name) {
    int idx = getTimerIdx(name);
    double localavg = times[idx]/counts[idx];
    return PCU_Max_Double(localavg);
  }

  double MigrationTimers::processAvg(std::string name) {
    int idx = getTimerIdx(name);
    double sum = PCU_Add_Double(times[idx]);
    return sum/PCU_Comm_Peers();
  }

  double MigrationTimers::perCallProcessAvg(std::string name) {
    int idx = getTimerIdx(name);
    double localavg = times[idx]/counts[idx];
    double sum = PCU_Add_Double(localavg);
    return sum/PCU_Comm_Peers();
  }

  int MigrationTimers::getTimerIdx(std::string name) {
    assert( nameToIdx.count(name) );
    return nameToIdx[name];
  }
} //end namespace
