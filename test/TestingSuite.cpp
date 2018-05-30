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

int TestingSuite::runTests() const {
  //Accumulation of failed tests
  int failures = 0;
  int num_tests = 0;
  //Run fine grain tests
  for (unsigned int i=0;i<fine_tests.size();i++) {
    int ierr = fine_tests[i]();
    num_tests++;
    if (ierr != 0) {
      EnGPar_Error_Message("Fine Test %d failed with error code: %d\n", i, ierr);
      failures++;
    }
  }

  //Run general tests on each graph
  for (unsigned int i = 0; i < general_tests.size(); i++) {
    for (unsigned int j = 0;j < test_graphs->size(); j++) {
      int ierr = general_tests[i](test_graphs->at(j));
      num_tests++;
      if (ierr != 0) {
        EnGPar_Error_Message("General Test %d failed with error code: %d\n", i, ierr);
        failures++;
      }
    }
  }

  if (!PCU_Comm_Self()) {
    EnGPar_Status_Message(-1,"%d/%ld tests passed in %s\n",num_tests-failures,
                          num_tests, n);
    if (failures==0) {
      EnGPar_Status_Message(-1,"All Tests Passed\n");
    }
  }
  return (failures==0);
}

