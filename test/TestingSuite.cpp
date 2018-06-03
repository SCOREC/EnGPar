#include "TestingSuite.h"
#include <engpar_support.h>
TestingSuite::TestingSuite(char* name) {
  n = name;
}
void TestingSuite::deleteTestGraphs() {
  for (unsigned int i=0;i<test_graphs.size();i++)
    agi::destroyGraph(test_graphs[i]);
  test_graphs.clear();
}

void TestingSuite::addFineTest(std::string name, fineTest t) {
  fine_names.push_back(name);
  fine_tests.push_back(t);
}
void TestingSuite::addGeneralTest(std::string name, generalTest t) {
  general_names.push_back(name);
  general_tests.push_back(t);
}
void TestingSuite::addTestGraph(std::string name, agi::Ngraph* g) {
  graph_names.push_back(name);
  test_graphs.push_back(g);
}

void TestingSuite::fillEmptyTestGraphs() {
  unsigned int max_graphs = PCU_Max_Int(test_graphs.size());
  for (unsigned int i=0;i<max_graphs;i++) {
    if (test_graphs.size()<=i) {
      test_graphs.push_back(agi::createEmptyGraph());
      graph_names.push_back("");
    }
  }
}

int TestingSuite::runTests(int trial) const {
  //Accumulation of failed tests
  int failures = 0;
  int num_tests = 0;
  //Run fine grain tests
  for (unsigned int i=0;i<fine_tests.size();i++) {
    int ierr = 0;
    int fail = 0;
    if (trial==-1||num_tests==trial) {
      if (!PCU_Comm_Self())
        EnGPar_Status_Message(-1,"Running Test %d: \"%s\".\n", num_tests,fine_names[i].c_str());
      EnGPar_Start_Test();
      ierr = fine_tests[i]();
      PCU_Barrier();      
    }
    if (ierr != 0) {
      EnGPar_Error_Message("Test %d: \"%s\"failed with error code: %d.\n", num_tests,
                           fine_names[i].c_str(), ierr);
      fail=1;
    }
    else if (EnGPar_Hit_Error()) {
      EnGPar_Error_Message("Test %d: \"%s\" failed with error message.\n", num_tests,
                           fine_names[i].c_str());
      fail = 1;
    }
    failures += PCU_Max_Int(fail);
    
    num_tests++;
  }

  //Run general tests on each graph
  for (unsigned int i = 0; i < general_tests.size(); i++) {
    for (unsigned int j = 0;j < test_graphs.size(); j++) {
      int ierr = 0;
      int fail = 0;
      if (trial==-1||num_tests==trial) {
        if (!PCU_Comm_Self())
          EnGPar_Status_Message(-1,"Running test %d: \"%s\" with graph: %s.\n", num_tests,
                                general_names[i].c_str(), graph_names[j].c_str());
        EnGPar_Start_Test();
        ierr = general_tests[i](test_graphs[j]);
        PCU_Barrier();
      }
      if (ierr != 0) {
        EnGPar_Error_Message("Test %d: \"%s\" with graph %s failed with error code: %d.\n",
                             num_tests, general_names[i].c_str(), graph_names[j].c_str(), ierr);
        fail = 1;
      }
      else if (EnGPar_Hit_Error()) {
        EnGPar_Error_Message("Test %d: \"%s\" failed with error message.\n", num_tests,
                             general_names[i].c_str());
        fail = 1;
      }
      failures += PCU_Max_Int(fail);
      num_tests++;
    }
  }
  PCU_Barrier();

  if (!PCU_Comm_Self()&&trial==-1) {
    EnGPar_Status_Message(-1,"%d/%ld tests passed in %s\n",num_tests-failures,
                          num_tests, n);
    if (failures==0) 
      EnGPar_Status_Message(-1,"All Tests Passed\n");
  }
  return failures;
}

