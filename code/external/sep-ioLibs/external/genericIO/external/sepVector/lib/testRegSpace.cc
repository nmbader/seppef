#include <gtest/gtest.h>  // googletest header file
#include <string>
#include "float2DReg.h"
#include "regSpace.h"
using std::string;
using namespace SEP;

TEST(regspace, axisToKey) {
  axis axisT(100, 0, 1.), axis2(10, 1., .5);
  std::shared_ptr<hypercube> hyper(new hypercube(axisT, axis2));
  float2DReg f(hyper);
  std::vector<float> k = f.axisToKey(2);

  ASSERT_EQ(k.size(), 10);
  for (int i = 0; i < 10; i++) ASSERT_EQ(k[i], 1 + .5 * i);
}