#include "TestingSuite.h"
#include <engpar_support.h>
TestingSuite::TestingSuite(char* name) {
  n = name;
}
void TestingSuite::deleteTestGraphs() {
  for (unsigned int i=0;i<test_graphs->size();i++)
    agi::destroyGraph(test_graphs->at(i));
  test_graphs->clear();
}

void TestingSuite::addFineTest(std::string name, fineTest t) {
  fine_names.push_back(name);
  fine_tests.push_back(t);
}
void TestingSuite::addGeneralTest(std::string name, generalTest t) {
  general_names.push_back(name);
  general_tests.push_back(t);
}
void TestingSuite::setTestingGraphs(std::vector<agi::Ngraph*>* gs) {
  test_graphs = gs;
}

int TestingSuite::runTests(int trial) const {
  //Accumulation of failed tests
  int failures = 0;
  int num_tests = 0;
  //Run fine grain tests
  for (unsigned int i=0;i<fine_tests.size();i++) {
    int ierr = 0;
    if (trial==-1||num_tests==trial) {
      EnGPar_Status_Message(-1,"Running Test %d: \"%s\"\n", num_tests,fine_names[i].c_str());
      ierr = fine_tests[i]();
      PCU_Barrier();
    }
    if (ierr != 0) {
      EnGPar_Error_Message("Test %d: \"%s\"failed with error code: %d\n", num_tests,
                           fine_names[i].c_str(), ierr);
      failures++;
    }
    num_tests++;
  }

  //Run general tests on each graph
  for (unsigned int i = 0; i < general_tests.size(); i++) {
    for (unsigned int j = 0;j < test_graphs->size(); j++) {
      int ierr = 0;
      if (trial==-1||num_tests==trial) {
        EnGPar_Status_Message(-1,"Running test %d: \"%s\" with graph: %d\n", num_tests,
                              general_names[i].c_str(), j);
        ierr = general_tests[i](test_graphs->at(j));
        PCU_Barrier();
      }
      if (ierr != 0) {
        EnGPar_Error_Message("Test %d: \"%s\" with graph %d failed with error code: %d\n",
                             num_tests, general_names[i].c_str(),j, ierr);
        failures++;
      }
      num_tests++;
    }
  }

  if (!PCU_Comm_Self()&&trial==-1) {
    EnGPar_Status_Message(-1,"%d/%ld tests passed in %s\n",num_tests-failures,
                          num_tests, n);
    if (failures==0) 
      EnGPar_Status_Message(-1,"All Tests Passed\n");
  }
  return failures;
}

