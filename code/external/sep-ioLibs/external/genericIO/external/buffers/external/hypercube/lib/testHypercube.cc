#include <gtest/gtest.h>  // googletest header file

#include <string>
#include "hypercube.h"
using std::string;
using namespace SEP;

void testHyper(std::shared_ptr<hypercube> h) {
  int nd = h->getNdim();
  for (int i = 0; i < nd; i++) {
    axis a = h->getAxis(i + 1);
    EXPECT_EQ(a.n, (i + 1) * 2);
    EXPECT_EQ(a.o, 3. * i);
    EXPECT_EQ(a.d, 4. * i);
    EXPECT_EQ(a.label, std::to_string(i));
    EXPECT_EQ(a.unit, std::to_string(i * 2));
  }
}

std::shared_ptr<hypercube> createHyper(const int ndim) {
  axis a1 = axis(2, 0., 0., std::to_string(0), std::to_string(0));
  if (ndim == 1) return std::make_shared<hypercube>(a1);
  axis a2 = axis(4, 3., 4., std::to_string(1), std::to_string(2));
  if (ndim == 2) return std::make_shared<hypercube>(a1, a2);
  axis a3 = axis(6, 6., 8., std::to_string(2), std::to_string(4));
  if (ndim == 3) return std::make_shared<hypercube>(a1, a2, a3);
  axis a4 = axis(8, 9., 12., std::to_string(3), std::to_string(6));
  if (ndim == 4) return std::make_shared<hypercube>(a1, a2, a3, a4);
  axis a5 = axis(10, 12., 16., std::to_string(4), std::to_string(8));
  if (ndim == 5) return std::make_shared<hypercube>(a1, a2, a3, a4, a5);
  axis a6 = axis(12, 15., 20., std::to_string(5), std::to_string(10));
  return std::make_shared<hypercube>(a1, a2, a3, a4, a5, a6);
}

TEST(hypercubeCreate, basic1D) { testHyper(createHyper(1)); }
TEST(hypercubeCreate, basic2D) { testHyper(createHyper(2)); }

TEST(hypercubeCreate, basic3D) { testHyper(createHyper(3)); }

TEST(hypercubeCreate, basic4D) { testHyper(createHyper(4)); }

TEST(hypercubeCreate, basic5D) { testHyper(createHyper(5)); }
TEST(hypercubeCreate, basic6D) { testHyper(createHyper(6)); }
TEST(hypercubeCreate, basicAxes) {
  std::shared_ptr<hypercube> hyp = createHyper(6),
                             hyp2 = std::make_shared<hypercube>(hyp->getAxes());
  testHyper(hyp2);
}
TEST(hypercubeCreate, clone) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  testHyper(hyp2);
}
TEST(hypercubeCheckSame, equal) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  ASSERT_NO_THROW(hyp2->checkSame(hyp));
}
TEST(hypercubeCheckSame, nwrong) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  axis a = hyp2->getAxis(4);
  a.n += 1;
  hyp2->setAxis(4, a);
  ASSERT_ANY_THROW(hyp2->checkSame(hyp));
}
TEST(hypercubeCheckSame, owrong) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  axis a = hyp2->getAxis(4);
  a.o += 1;
  hyp2->setAxis(4, a);
  ASSERT_ANY_THROW(hyp2->checkSame(hyp));
}
TEST(hypercubeCheckSame, dwrong) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  axis a = hyp2->getAxis(4);
  a.d += 1;
  hyp2->setAxis(4, a);
  ASSERT_ANY_THROW(hyp2->checkSame(hyp));
}
TEST(hypercubeCheckSame, olittlewrong) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  axis a = hyp2->getAxis(4);
  a.o += .0001;
  hyp2->setAxis(4, a);
  ASSERT_NO_THROW(hyp2->checkSame(hyp));
}

TEST(hypercubeCheckSame, dlittlewrong) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  axis a = hyp2->getAxis(4);
  a.d += .0001;
  hyp2->setAxis(4, a);
  ASSERT_NO_THROW(hyp2->checkSame(hyp));
}

TEST(hypercubeCheckSameSize, equal) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  EXPECT_TRUE(hyp2->sameSize(hyp));
}
TEST(hypercubeCheckSameSize, nwrong) {
  std::shared_ptr<hypercube> hyp = createHyper(6), hyp2 = hyp->clone();
  axis a = hyp2->getAxis(4);
  a.n += 1;
  hyp2->setAxis(4, a);
  EXPECT_FALSE(hyp2->sameSize(hyp));
}
TEST(hypercubeSizes, get123) {
  std::shared_ptr<hypercube> hyp = createHyper(6);
  EXPECT_EQ(hyp->getN123(), 46080);
}
TEST(hypercubeSizes, getNs) {
  std::shared_ptr<hypercube> hyp = createHyper(6);
  for (int i = 0; i < 6; i++) {
    EXPECT_EQ(hyp->getAxis(i + 1).n, (i + 1) * 2);
  }
}
TEST(hypercubeCreate, setAxes) {
  std::shared_ptr<hypercube> hyp = createHyper(1), hyp2 = createHyper(6);
  hyp->setAxes(hyp2->getAxes());
  testHyper(hyp);
}

TEST(hypercubeSizes, ndim) {
  std::shared_ptr<hypercube> hyp = createHyper(6);
  ASSERT_EQ(6, hyp->getNdim());
  axis a = axis(1);
  hyp->addAxis(a);
  ASSERT_EQ(6, hyp->getNdimG1());
  ASSERT_EQ(7, hyp->getNdim());
}
TEST(hypercubeSizes, shrink) {
  std::shared_ptr<hypercube> hyp = createHyper(6);
  ASSERT_EQ(6, hyp->getNdim());
  axis a = axis(4);

  hyp->addAxis(a);
  ASSERT_EQ(7, hyp->getNdim());
  hyp->shrinkDimension(6);
  ASSERT_EQ(6, hyp->getNdim());
}
