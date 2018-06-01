#ifndef TESTING_SUITE_H__
#define TESTING_SUITE_H__
#include <vector>
#include <ngraph.h>
class TestingSuite {
 public:
  //Fine-grain test functions
  //  returns an error code (0 = success)
  typedef int (*fineTest)();
  //General test functions
  //  returns an error code (0 = success)
  typedef int (*generalTest)(agi::Ngraph*);

  TestingSuite(char* name);
  void deleteTestGraphs();
  
  void addFineTest(fineTest t);
  void addGeneralTest(generalTest t);
  void setTestingGraphs(std::vector<agi::Ngraph*>* gs);

  int runTests(int trial = -1) const;

 private:
  char* n;
  std::vector<fineTest> fine_tests;
  std::vector<generalTest> general_tests;
  std::vector<agi::Ngraph*>* test_graphs;
};
#endif
