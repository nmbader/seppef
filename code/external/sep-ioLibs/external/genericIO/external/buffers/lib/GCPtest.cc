#include <gtest/gtest.h>  // googletest header file
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include "gcpBuffers.h"
#include "google/cloud/status_or.h"
#include "google/cloud/storage/client.h"
#include "google/cloud/storage/oauth2/google_credentials.h"
#include "ioTypes.h"
#include "nocompress.h"
using namespace std::chrono;
using ::google::cloud::StatusOr;

using std::string;
using namespace SEP::IO;

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

/*
TEST(TESTBucketCreation, smallTest) {

  long long n123 = 1;
  int ndim = 3;
  int n=80;
  std::vector<int> ns(ndim, n), fs(ndim, 0), js(ndim, 1);


  std::string bucket = std::string("test-junk-1");
  std::string test="test";
  std::shared_ptr<noCompression> comp(new noCompression(SEP::DATA_FLOAT));

  // Create simple file and write to disk
  gcpBuffer buf(bucket,ns,fs,comp,UNDEFINED);
  buf.setName(test);
  buf.changeState(CPU_DECOMPRESSED);
  std::shared_ptr<storeFloat> inA=array();

  buf.putBufferCPU(inA,CPU_DECOMPRESSED);
  buf.changeState(ON_DISK);
  std::shared_ptr<storeBase> outA=inA->clone();
  sleep(5);
  buf.changeState(CPU_DECOMPRESSED);
  buf.getBufferCPU(outA,CPU_DECOMPRESSED);

  float *inp=(float*)inA->getPtr(),*outp=(float*)outA->getPtr();
  for(int i=0;i <27; i++){
           EXPECT_EQ(inp[i*13],outp[i*13]);
  }


}
*/
TEST(TESTGCP, basic) {
  namespace gcs = google::cloud::storage;

  google::cloud::v0::StatusOr<gcs::Client> client2 =
      gcs::Client::CreateDefaultClient();

  EXPECT_TRUE(client2);
  long desired_line_count = 10;
  std::string bucket_name = "unit-test-1234";
  std::string object_name = "SSSS";

  std::string projectID =
      SEP::IO::getEnvVar(std::string("projectID"), std::string("earth-clapp"));
  std::string region =
      SEP::IO::getEnvVar(std::string("region"), std::string("us-west1"));
  client2.value().DeleteBucket(bucket_name);

  google::cloud::StatusOr<gcs::BucketMetadata> bucket_metadata =
      client2->CreateBucketForProject(
          bucket_name, projectID,
          gcs::BucketMetadata().set_location(region).set_storage_class(
              gcs::storage_class::Regional()));

  EXPECT_TRUE(bucket_metadata);
  std::cerr << bucket_metadata.status().message() << std::string(":")
            << std::endl;
  std::string const text = "Lorem ipsum dolor sit amet";
  gcs::ObjectWriteStream stream =
      client2->WriteObject(bucket_name, object_name);
  for (int lineno = 0; lineno != desired_line_count; ++lineno) {
    // Add 1 to the counter, because it is conventional to number lines
    // starting at 1.
    stream << (lineno + 1) << ": " << text << "\n";
  }
  stream.Close();
  StatusOr<gcs::ObjectMetadata> metadata = std::move(stream).metadata();
  EXPECT_TRUE(metadata);
  client2.value().DeleteObject(bucket_name, object_name);

  google::cloud::Status status = client2.value().DeleteBucket(bucket_name);

  if (!status.ok()) {
    throw std::runtime_error(status.message());
  }

  // client2.value().DeleteBucket(bucket_name);
}

TEST(TESTGCP, basicBuffer) {
  std::vector<SEP::axis> axes;
  long long n = 200;
  long long n123 = 1;
  int ndim = 4;
  std::vector<int> ns(ndim, n), fs(ndim, 0), js(ndim, 1);
  for (int i = 0; i < ndim; i++) {
    n123 = n123 * n;
    axes.push_back(SEP::axis(n));
  }

  std::shared_ptr<SEP::hypercube> hyper(new SEP::hypercube(axes));

  std::string bucket = std::string("testbucket99r");
  std::string bucket1 = bucket + std::string("/dataset1");
  std::string bucket2 = bucket + std::string("/dataset2");

  // Create Default blocking
  //  std::shared_ptr<SEP::IO::blocking> block =
  //      SEP::IO::blocking::createDefaultBlocking(hyper);

  std::vector<int> big(4, 40), bs(4, 2);
  big[0] = 100;
  std::shared_ptr<SEP::IO::blocking> block(new SEP::IO::blocking(bs, big));

  float *vals = new float[n123];

  Json::Value val;
  high_resolution_clock::time_point t2, t3, t1;
  int ntimes = 1;
  for (int i = 0; i < ntimes; i++) {
    // Create simple file and write to disk
    SEP::IO::gcpBuffers gcp(hyper, SEP::DATA_FLOAT, block);
    ASSERT_NO_THROW(gcp.setName(bucket1 + std::to_string(i), true));

    //  std::vector<float> vals(n123);

    t1 = high_resolution_clock::now();

    ASSERT_NO_THROW(gcp.putWindow(ns, fs, js, vals));
    t2 = high_resolution_clock::now();

    ASSERT_NO_THROW(gcp.changeState(ON_DISK));
    high_resolution_clock::time_point t3 = high_resolution_clock::now();

    auto d1 = duration_cast<microseconds>(t2 - t1).count();
    auto d2 = duration_cast<microseconds>(t3 - t2).count();
    double s1 = (double)n123 * 4 / d1;
    double s2 = (double)n123 * 4 / d2;

    std::cerr << "To buffer " << s1 << " MB/s " << std::endl;
    ;
    std::cerr << "To cloud " << s2 << " MB/s " << std::endl;
    val["blocking"] = block->getJsonDescription();
    val["compression"] = gcp.getCompressObj()->getJsonDescription();
  }

  // Create a second directory in same bucket

  // SEP::IO::gcpBuffers gcp2(hyper, SEP::DATA_FLOAT, block);
  // ASSERT_NO_THROW(gcp2.setName(bucket2, true));

  //  ASSERT_NO_THROW(gcp2.putWindow(ns, fs, js, vals.data()));
  //  ASSERT_NO_THROW(gcp2.changeState(ON_DISK));

  // Now read the bucket from disk

  float tot = 0;
  for (int i = 0; i < ntimes; i++) {
    float *vals2 = new float[n123];
    SEP::IO::gcpBuffers gcp3(hyper, bucket1 + std::to_string(i), val);
    t2 = high_resolution_clock::now();

    //  ASSERT_NO_THROW(gcp3.getWindow(ns, fs, js, vals2.data()));
    //  ASSERT_NO_THROW(gcp3.changeState(CPU_COMPRESSED));
    ASSERT_NO_THROW(gcp3.getWindow(ns, fs, js, vals2));

    t3 = high_resolution_clock::now();
    auto d2 = duration_cast<microseconds>(t3 - t2).count();
    std::cerr << "From cloud " << (double)n123 * 4 / d2 << " MB/s "
              << std::endl;
    ;
    tot += d2;
  }

  std::cerr << 4 * n123 * ntimes / tot << " average speed" << std::endl;

  namespace gcs = google::cloud::storage;
  google::cloud::v0::StatusOr<gcs::Client> client =
      gcs::Client::CreateDefaultClient();

  for (auto &&object_metadata : client.value().ListObjects(bucket)) {
    if (!object_metadata) {
      throw std::runtime_error(object_metadata.status().message());
    }
    client.value().DeleteObject(bucket, object_metadata->name());
  }
  client.value().DeleteBucket(bucket);
}
