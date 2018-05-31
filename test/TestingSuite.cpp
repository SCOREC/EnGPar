#include "TestingSuite.h"
#include <engpar_support.h>
TestingSuite::TestingSuite(char* name) {
  n = name;
}

void TestingSuite::addFineTest(fineTest t) {
  fine_tests.push_back(t);
}
void TestingSuite::addGeneralTest(generalTest t) {
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
    if (trial==-1||num_tests==trial)
      ierr = fine_tests[i]();
    if (ierr != 0) {
      EnGPar_Error_Message("Test: %d, Fine Test %d failed with error code: %d\n", num_tests,i, ierr);
      failures++;
    }
    num_tests++;
  }

  //Run general tests on each graph
  for (unsigned int i = 0; i < general_tests.size(); i++) {
    for (unsigned int j = 0;j < test_graphs->size(); j++) {
      int ierr = 0;
      if (trial==-1||num_tests==trial)
        ierr = general_tests[i](test_graphs->at(j));
      if (ierr != 0) {
        EnGPar_Error_Message("Test: %d, General Test %d Graph %d failed with error code: %d\n", num_tests, i, j, ierr);
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
  return (failures==0);
}

