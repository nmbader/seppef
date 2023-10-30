#include "gcpBuffersRegFile.h"
#include "gcpBuffers.h"
#include "google/cloud/storage/client.h"
#include "google/cloud/storage/oauth2/google_credentials.h"
#include <cstdlib>
#include <exception>
#include <fstream>  // std::ifstream
#include <iostream> // std::cout
using namespace SEP;
gcpBuffersRegFile::gcpBuffersRegFile(const Json::Value &arg,
                                     const usage_code usage,
                                     const std::string &tag,
                                     const std::string &progName,
                                     const int ndimMax) {
  setUsage(usage);
  setupGCP(arg, tag);
  if (!_newFile) {
    readDescription(ndimMax);

    if (jsonArgs["bufferInfo"].isNull())
      error(std::string("bufferInfo not provided in JSON file"));
    _bufs.reset(
        new SEP::IO::gcpBuffers(getHyper(), tag, jsonArgs["bufferInfo"]));
  }
  _usage = usage;
  jsonArgs["progName"] = progName;
  jsonArgs["name"] = tag;
  _type = "GCP Buffers";
  _binary = "Multiple objects";
}

void gcpBuffersRegFile::setupGCP(const Json::Value &arg,
                                 const std::string &tag) {
  _tag = tag;
  std::string bucket, baseName;

  int pos;
  if ((pos = tag.find("/")) == std::string::npos) { // No subdirectory
    bucket = tag;
    baseName = "";
  } else {
    baseName = tag;
    bucket = baseName.substr(0, baseName.find("/"));
    baseName.erase(0, baseName.find("/") + 1);
  }
  _bucket = bucket;
  _dir = baseName;
  std::string _projectID = getEnvVar("projectID", "NONE");
  std::string _region = getEnvVar("region", "us-west1");
  if (_projectID == std::string("NONE")) {
    std::cerr << "Must set environmental variable projectID " << _projectID
              << std::endl;
    exit(1);
  }
  namespace gcs = google::cloud::storage;
  google::cloud::v0::StatusOr<gcs::Client> client =
      gcs::Client::CreateDefaultClient();
  if (!client)
    throw(SEPException(std::string("Trouble creating default client")));
  _client = client;
  if (arg[tag].isNull()) {
    _jsonFile = _tag;
  } else {
    _jsonFile = arg[tag].asString();
  }
  _newFile = true;
  if (_usage == usageIn)
    _newFile = false;
  else if (_usage == usageInOut) {
    try {
      [&](gcs::Client client, std::string bucket_name,
          std::string object_name) {
        gcs::ObjectReadStream stream =
            client.ReadObject(bucket_name, object_name);
        std::string data(std::istreambuf_iterator<char>(stream), {});
        std::string errs;
        Json::CharReaderBuilder builder;
        Json::parseFromStream(builder, stream, &jsonArgs, &errs);
        //        jsonArgs=JSON::parse(data);
        stream.Close();
      }(std::move(_client.value()), _bucket, _dir + std::string("/desc"));
    } catch (std::exception const &ex) {
      std::cerr << "Trouble reading from bucket " << _bucket + _dir + "/desc"
                << std::endl;
      exit(1);
    }
  }
  if (_usage == usageIn || !_newFile) {
    try {
      [&](gcs::Client client, std::string bucket_name,
          std::string object_name) {
        gcs::ObjectReadStream stream =
            client.ReadObject(bucket_name, object_name);
        std::string data(std::istreambuf_iterator<char>(stream), {});
        Json::Reader read;
        read.parse(data, jsonArgs);
        stream.Close();
      }(std::move(_client.value()), _bucket, _dir + std::string("/desc"));
    } catch (std::exception const &ex) {
      std::cerr << "Trouble reading from bucket "
                << _bucket + _dir + std::string("/desc") << std::endl;
      exit(1);
    }
  }
}
void gcpBuffersRegFile::removeDescDir() {
  namespace gcs = google::cloud::storage;
  google::cloud::v0::StatusOr<gcs::Client> client =
      gcs::Client::CreateDefaultClient();

  client.value().DeleteObject(_bucket, _dir + std::string("/desc"));
}
void gcpBuffersRegFile::close() {
  writeDescription();
  namespace gcs = google::cloud::storage;
  if (getUsage() == usageOut || getUsage() == usageInOut) {
    gcs::ObjectWriteStream stream =
        _client.value().WriteObject(_bucket, _dir + std::string("/desc"));
    stream << jsonArgs;
    stream.Close();
    google::cloud::v0::StatusOr<gcs::ObjectMetadata> metadata =
        std::move(stream).metadata();
    if (!metadata) {
      std::cerr << "trouble closing " << _bucket + std::string("/") + _dir
                << std::endl;
      std::cerr << metadata.status().message() << std::endl;
      throw SEPException(std::string("Trouble writing object"));
    }
  }
  _bufs->changeState(SEP::IO::ON_DISK);
}
void gcpBuffersRegFile::createBuffers() {
  if (_bufs) {
    return;
  }
  if (!_hyper)
    error("Must set hypercube before blocking");
  if (getDataType() == SEP::DATA_UNKNOWN)
    error("Must set dataType before setting blocks");

  _bufs.reset(
      new SEP::IO::gcpBuffers(getHyper(), getDataType(), _block, _comp, _mem));

  _bufs->setName(jsonArgs["name"].asString(), true);
}
