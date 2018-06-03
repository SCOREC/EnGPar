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
  
  void addFineTest(std::string name, fineTest t);
  void addGeneralTest(std::string name, generalTest t);
  void addTestGraph(std::string name, agi::Ngraph* g);
  void fillEmptyTestGraphs();
  
  int runTests(int trial = -1) const;

 private:
  char* n;
  std::vector<std::string> fine_names;
  std::vector<fineTest> fine_tests;
  std::vector<std::string> general_names;
  std::vector<generalTest> general_tests;
  std::vector<std::string> graph_names;
  std::vector<agi::Ngraph*> test_graphs;
  
};
#endif
