
#include <gtest/gtest.h>  // googletest header file
#include <fstream>
#include <iostream>
#include <string>
#include "ZfpCompress.h"
#include "fileBuffer.h"
#include "nocompress.h"
using std::string;
using namespace SEP::IO;
using namespace SEP;

TEST(writeBuffer, buffer) {
  std::shared_ptr<noCompression> comp(new noCompression(DATA_INT));

  std::vector<int> n(3, 20), f(3, 0);

  std::vector<int> a(20 * 20 * 20, 0);
  for (auto i = 0; i < 20 * 20 * 20; i++) a[i] = i;
  std::shared_ptr<storeInt> store(new storeInt(20 * 20 * 20, (void *)a.data()));

  fileBuffer buf(n, f, comp, UNDEFINED);

  std::vector<int> nblock = buf.getBlock();

  EXPECT_EQ(nblock[0], 1);
  EXPECT_EQ(nblock[1], 20);
  EXPECT_EQ(nblock[2], 400);

  buf.setName("/tmp/junk1");

  buf.putBufferCPU(store, CPU_DECOMPRESSED);

  ASSERT_NO_THROW(buf.writeBuffer());
  std::ifstream in("/tmp/junk1", std::ifstream::ate | std::ifstream::binary);

  ASSERT_NO_THROW(in.good());
  size_t nelemFile = in.tellg();
  in.seekg(0, std::iostream::beg);

  EXPECT_EQ(nelemFile, 20 * 20 * 20 * 4);
  std::vector<int> b(8000, 0.);
  in.read((char *)b.data(), nelemFile);

  in.close();
  for (int i = 0; i < 20 * 20 * 20; i += 29) {
    EXPECT_EQ(i, b[i]);
  }
}

std::shared_ptr<storeFloat> array() {
  std::vector<float> ar(8000);
  size_t ii = 0;

  for (auto i3 = 0; i3 < 20; i3++) {
    for (auto i2 = 0; i2 < 20; i2++) {
      for (auto i1 = 0; i1 < 20; i1++, ii++) {
        ar[ii] = (float)rand() / ((float)RAND_MAX) + i1 * .4 + 5.;
      }
    }
  }
  std::shared_ptr<storeFloat> store(new storeFloat(8000, ar.data()));
  return store;
}
TEST(readWindowNoCompress, buffer) {
  ZfpParams zpars = ZfpParams();
  zpars._meth = ZFP_RATE;
  std::shared_ptr<noCompression> comp(new noCompression(DATA_FLOAT));

  std::shared_ptr<storeFloat> store = array();
  std::vector<int> n(3, 20), f(3, 0);
  fileBuffer buf(n, f, comp, UNDEFINED);
  buf.putBufferCPU(store, CPU_DECOMPRESSED);
  buf.setName("/tmp/junk2");
  std::shared_ptr<storeBase> storeCompare = store->clone();

  ASSERT_NO_THROW(buf.writeBuffer());
  std::cerr<<buf._bufferState<<std::endl;
  buf.readBuffer();
  ASSERT_NO_THROW(buf.readBuffer());

  std::shared_ptr<storeFloat> x(new storeFloat(20 * 20 * 20));

  float *b1 = (float *)storeCompare->getPtr();
  buf.getBufferCPU(x, CPU_DECOMPRESSED);
  float *b2 = (float *)x->getPtr();

  for (int ii = 0; ii < 20 * 20 * 20; ii++) {
    ASSERT_LE(fabsf((b1[ii] - b2[ii]) / b2[ii]), .01);
  }
}

TEST(readWindowCompress, buffer) {
  ZfpParams zpars = ZfpParams();
  zpars._meth = ZFP_RATE;
  std::shared_ptr<ZfpCompression> z(new ZfpCompression(DATA_FLOAT, zpars));

  std::shared_ptr<storeFloat> store = array();
  std::vector<int> n(3, 20), f(3, 0);
  fileBuffer buf(n, f, z, UNDEFINED);
  buf.putBufferCPU(store, CPU_DECOMPRESSED);
  buf.setName("/tmp/junk2");
  std::shared_ptr<storeBase> storeCompare = store->clone();

  ASSERT_NO_THROW(buf.writeBuffer());
  ASSERT_NO_THROW(buf.readBuffer());

  std::shared_ptr<storeFloat> x(new storeFloat(20 * 20 * 20));

  float *b1 = (float *)storeCompare->getPtr();
  buf.getBufferCPU(x, CPU_DECOMPRESSED);
  float *b2 = (float *)x->getPtr();

  for (int ii = 0; ii < 20 * 20 * 20; ii++) {
    ASSERT_LE(fabsf((b1[ii] - b2[ii]) / b2[ii]), .01);
  }
}

/*


  zpars = ZfpParams();
  zpars._meth = ZFP_TOLERANCE;
  compare(zpars, .012, 3., 7.);
  */

void createAxisTest(const int _f, const int _n, const int nw, const int fw,
                    const int jw, const int f_w, const int fG, const int n) {
  std::shared_ptr<noCompression> comp(new noCompression(DATA_INT));

  std::vector<int> nbuf(7, _n), fbuf(7, _f);
  fileBuffer buf(nbuf, fbuf, comp, UNDEFINED);
  std::vector<int> nwS(1, nw), fwS(1, fw), jwS(1, jw), f_wR(1), fGR(1), nR(1),
      blockG(2), jl(1), ng(1);
  buf.localWindow(nwS, fwS, jwS, nR, f_wR, jl, ng, fGR, blockG);
  EXPECT_EQ(nR[0], n);
  EXPECT_EQ(f_wR[0], f_w);
  EXPECT_EQ(fGR[0], fG);
}

TEST(localWindow, buffer) {
  std::shared_ptr<noCompression> comp(new noCompression(DATA_INT));

  //              _f,  _n,   nw,  fw,   jw,   f_w,   fG,   n
  std::cerr << "test1" << std::endl;
  createAxisTest(11, 20, 100, 0, 1, 0, 11, 20);
  std::cerr << "test2" << std::endl;

  createAxisTest(11, 20, 100, 0, 3, 1, 4, 7);
  std::cerr << "test3" << std::endl;

  createAxisTest(240, 60, 248, 248, 1, 8, 0, 52);
}

/*
  buffer(const std::string name, const std::vector<int> &n,
         const std::vector<int> &f,
         std::shared_ptr<compress> comp);  // Read from file
  buffer(const std::vector<int> &n, const std::vector<int> &f,
         std::shared_ptr<compress> comp);


  long long getWindowCPU(const std::vector<int> &nw,
                         const std ::vector<int> &fw,
                         const std::vector<int> &jw,
                         std::shared_ptr<storeBase> buf, const size_t bufLoc,
                         const bool keepState = false);

  long long putWindowCPU(const std::vector<int> &nw,
                         const std ::vector<int> &fw, std::vector<int> &jw,
                         const std::shared_ptr<storeBase> buf,
                         const size_t bufLoc, const bool keepState = false);

  long changeState(const bufferState state);
  */
