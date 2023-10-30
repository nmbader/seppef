#include <gtest/gtest.h>  // googletest header file

#include <string>
#include "genericFile.h"
#include "ioModes.h"
#include "memoryFile.h"
using std::string;
using namespace SEP;
std::vector<float> createArrayF(const int n1, const int n2, const int n3) {
  long long n123 = (long long)n1 * (long long)n2 * (long long)n3;
  std::vector<float> buf(n123);
  long long i = 0;
  for (int i3 = 0; i3 < n3; i3++) {
    for (int i2 = 0; i2 < n2; i2++) {
      for (int i1 = 0; i1 < n1; i1++, i++) {
        buf[i] = i1 + i2 * 100 + i3 * 100 * 100;
      }
    }
  }
  return buf;
}
void checkArrayF(const float *buf, const int n1, const int f1, const int j1,
                 const int n2, const int f2, const int j2, const int n3,
                 const int f3, const int j3) {
  long long i = 0;
  for (int i3 = 0; i3 < n3; i3++) {
    for (int i2 = 0; i2 < n2; i2++) {
      for (int i1 = 0; i1 < n1; i1++, i++) {
        EXPECT_EQ(buf[i], (f1 + i1 * j1) + (f2 + i2 * j2) * 100 +
                              (f3 + j3 * i3) * 100 * 100);
      }
    }
  }
}

TEST(initial, getIO) {
  std::vector<std::string> args;
  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;

  ASSERT_NO_THROW(io = modes.getIO("memory"));
}

TEST(initial, createBufferf) {
  std::vector<std::string> args;
  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;

  io = modes.getIO("memory");

  std::shared_ptr<genericRegFile> file;
  ASSERT_NO_THROW(file = io->getRegFile("outputTest", usageOut));
}
TEST(initial, correctType) {
  std::vector<std::string> args;
  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;

  io = modes.getIO("memory");

  std::shared_ptr<genericRegFile> file = io->getRegFile("outputTest", usageOut);
  std::shared_ptr<memoryRegFile> flM =
      std::dynamic_pointer_cast<memoryRegFile>(file);
  EXPECT_NE(flM, nullptr);
}

TEST(writeFloat, hyperMustBeSet) {
  std::vector<std::string> args;
  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;

  io = modes.getIO("memory");

  std::shared_ptr<genericRegFile> file = io->getRegFile("outputTest", usageOut);
  std::vector<float> buf1 = createArrayF(10, 10, 10);

  ASSERT_ANY_THROW(file->writeFloatStream(buf1.data(), 1000));
}

TEST(writeFloat, All) {
  std::vector<std::string> args;

  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;
  io = modes.getIO("memory");
  std::shared_ptr<hypercube> hyper(new hypercube(10, 10, 10));

  std::shared_ptr<genericRegFile> file = io->getRegFile("outputTest", usageOut);
  file->setHyper(hyper);
  std::shared_ptr<memoryRegFile> flM =
      std::dynamic_pointer_cast<memoryRegFile>(file);

  std::vector<float> buf1 = createArrayF(10, 10, 10);

  file->writeFloatStream(buf1.data(), 1000);

  std::vector<unsigned char> v = flM->returnVec();

  EXPECT_EQ(v.size(), 4000);

  float *buf = (float *)flM->getPtr();

  checkArrayF(buf, 10, 0, 1, 10, 0, 1, 10, 0, 1);
}

TEST(writeFloat, parts) {
  std::vector<std::string> args;

  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;
  io = modes.getIO("memory");
  std::shared_ptr<hypercube> hyper(new hypercube(10, 10, 10));

  std::shared_ptr<genericRegFile> file = io->getRegFile("outputTest", usageOut);
  file->setHyper(hyper);
  std::shared_ptr<memoryRegFile> flM =
      std::dynamic_pointer_cast<memoryRegFile>(file);

  std::vector<float> buf1 = createArrayF(10, 10, 10);
  for (int i = 0; i < 10; i++)
    file->writeFloatStream(buf1.data() + 100 * i, 100);

  std::vector<unsigned char> v = flM->returnVec();

  EXPECT_EQ(v.size(), 4000);

  float *buf = (float *)flM->getPtr();

  checkArrayF(buf, 10, 0, 1, 10, 0, 1, 10, 0, 1);
}

TEST(writeFloat, wind) {
  std::vector<std::string> args;

  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;
  io = modes.getIO("memory");
  std::shared_ptr<hypercube> hyper(new hypercube(10, 10, 10));

  std::shared_ptr<genericRegFile> file = io->getRegFile("outputTest", usageOut);
  file->setHyper(hyper);
  std::shared_ptr<memoryRegFile> flM =
      std::dynamic_pointer_cast<memoryRegFile>(file);

  std::vector<float> buf1 = createArrayF(10, 10, 10);
  std::vector<int> nw(3, 10), fw(3, 0), jw(3, 1);
  nw[0] = 2;
  jw[0] = 5;
  std::vector<float> buf2(100 * 2);
  for (int i = 0; i < 5; i++) {
    fw[0] = i;
    for (int ii = 0; ii < 100; ii++) buf2[ii * 2] = buf1[ii * 10 + i];
    for (int ii = 0; ii < 100; ii++) buf2[ii * 2 + 1] = buf1[ii * 10 + i + 5];
    file->writeFloatWindow(nw, fw, jw, buf2.data());
  }
  float *buf = (float *)flM->getPtr();

  std::vector<unsigned char> v = flM->returnVec();

  EXPECT_EQ(v.size(), 4000);

  checkArrayF(buf, 10, 0, 1, 10, 0, 1, 10, 0, 1);
}

TEST(readFloat, All) {
  std::vector<std::string> args;

  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;
  io = modes.getIO("memory");
  std::shared_ptr<hypercube> hyper(new hypercube(10, 10, 10));

  std::shared_ptr<genericRegFile> file = io->getRegFile("outputTest", usageOut);
  file->setHyper(hyper);
  std::shared_ptr<memoryRegFile> flM =
      std::dynamic_pointer_cast<memoryRegFile>(file);

  std::vector<float> buf1 = createArrayF(10, 10, 10);

  file->writeFloatStream(buf1.data(), 1000);

  std::vector<float> buf2(1000, 0.);
  file->seekTo(0, 0);
  file->readFloatStream(buf2.data(), 1000);

  checkArrayF(buf2.data(), 10, 0, 1, 10, 0, 1, 10, 0, 1);
}

TEST(readFloat, Window) {
  std::vector<std::string> args;

  ioModes modes = ioModes(args);
  std::shared_ptr<genericIO> io;
  io = modes.getIO("memory");
  std::shared_ptr<hypercube> hyper(new hypercube(10, 10, 10));

  std::shared_ptr<genericRegFile> file = io->getRegFile("outputTest", usageOut);
  file->setHyper(hyper);
  std::shared_ptr<memoryRegFile> flM =
      std::dynamic_pointer_cast<memoryRegFile>(file);

  std::vector<float> buf1 = createArrayF(10, 10, 10);

  file->writeFloatStream(buf1.data(), 1000);

  std::vector<int> nw(3, 10), fw(3, 0), jw(3, 1);
  nw[0] = 2;
  jw[0] = 5;
  std::vector<float> buf2(100 * 2);
  for (int i = 0; i < 5; i++) {
    fw[0] = i;
    file->readFloatWindow(nw, fw, jw, buf2.data());

    for (int ii = 0; ii < 100; ii++) EXPECT_EQ(buf2[ii * 2], buf1[ii * 10 + i]);
    for (int ii = 0; ii < 100; ii++)
      EXPECT_EQ(buf2[ii * 2 + 1], buf1[ii * 10 + i + 5]);
  }
}