#include <gtest/gtest.h>  // googletest header file

#include <string>
#include "dictParams.h"
using std::string;
using namespace SEP;
std::map<std::string, std::string> createPars() {
  std::map<std::string, std::string> ps;
  ps["testI1"] = std::to_string(1);
  ps["testf1"] = std::to_string(2.);
  ps["tests1"] = "blah";
//  ps["testl1"] = false;
  ps["testiA"] = "2,3,4,6,8";

  ps["testfA"] = "2.5,3.5,4.5,5.5";
  return ps;
}

TEST(dictParams, intTest) {
  dictParams pars = dictParams(createPars());
  EXPECT_EQ(pars.getInt("testI1"), 1);
}
TEST(dictParams, intTest2) {
  dictParams pars = dictParams(createPars());
  EXPECT_EQ(pars.getInt("testI2", 2), 2);
}
TEST(dictParams, intTest3) {
  dictParams pars = dictParams(createPars());
  ASSERT_ANY_THROW(pars.getInt("testI2"));
}
TEST(dictParams, intTest4) {
  dictParams pars = dictParams(createPars());
  EXPECT_EQ(pars.getInt("testI1", 2), 1);
}
