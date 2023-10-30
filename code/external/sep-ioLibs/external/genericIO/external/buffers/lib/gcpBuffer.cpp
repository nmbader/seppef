#include <gcpBuffer.h>
#include <locale.h>
#include <stdio.h>
#include <unistd.h>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "google/cloud/status_or.h"
#include "google/cloud/storage/client.h"
#include "google/cloud/storage/oauth2/google_credentials.h"
//#include "statusor.h"
using namespace SEP::IO;

long long gcpBuffer::writeBuffer(bool keepState) {
  long long oldSize = _buf->getSize();

  std::shared_ptr<storeBase> buf;
  bufferState restore;
  if (_bufferState == UNDEFINED) throw SEPException("Bufferstate is undefined");
  if (_bufferState == ON_DISK) return 0;
  if (keepState) {
    restore = _bufferState;
    buf = _buf->clone();
  }

  changeState(CPU_COMPRESSED);
  namespace gcs = google::cloud::storage;
  //  google::cloud::v0::StatusOr<gcs::Client> client =
  //     gcs::Client::CreateDefaultClient();
  google::cloud::v0::StatusOr<gcs::ObjectMetadata> metadata;
  gcs::ObjectWriteStream stream;
  int sleep = 100000;  // Start with a retry at .1 seconds
  bool success = false;
  int itry = 0;

  while (!success && itry < _ntrys) {
    stream = _client.value().WriteObject(_bucketName, _name);

    stream.write(_buf->getPtr(), _buf->getSize() * _buf->getElementSize());

    stream.Close();

    metadata = std::move(stream).metadata();

    if (metadata) {
      success = true;
    } else {
      usleep(sleep);
      sleep = sleep * 3;
    }
    itry++;
  }

  if (!metadata) {
    std::cerr << "FAIL: " << _name << "after " << itry
              << " trys message=" << metadata.status().message()
              << std::string(":") << std::endl;
    throw SEPException(std::string("Trouble writing object"));
  }
  if (stream.received_hash() != stream.computed_hash())
    throw SEPException(std::string("Hashes do not match"));

  if (keepState) {
    _buf = buf;
    _bufferState = restore;
  } else {
    _bufferState = ON_DISK;
  }
  return _buf->getSize() - oldSize;
}

long long gcpBuffer::readBuffer() {
  long long oldSize = _buf->getSize();
  _modified = false;
  namespace gcs = google::cloud::storage;
  google::cloud::v0::StatusOr<gcs::ObjectMetadata> object_metadata;
  // google::cloud::v0::StatusOr<gcs::Client> client =
  //   gcs::Client::CreateDefaultClient();
  /*Only need to do something if sitting on disk*/

  namespace gcs = google::cloud::storage;
  if (_bufferState == ON_DISK) {
    std::shared_ptr<storeByte> buf = std::dynamic_pointer_cast<storeByte>(_buf);

    int sleep = 100000;  // Start with a retry at .1 seconds
    bool success = false;
    int itry = 0;
    while (!success && itry < _ntrys) {
      object_metadata = _client.value().GetObjectMetadata(_bucketName, _name);

      if (object_metadata)
        success = true;
      else {
        usleep(sleep);
        sleep = sleep * 3;
      }
      itry++;
    }
    if (!object_metadata) {
      throw std::runtime_error(object_metadata.status().message());
    }

    auto sz = object_metadata->size();

    gcs::ObjectReadStream stream =
        _client.value().ReadObject(_bucketName, _name);

    // if (!stream.IsOpen())
    // throw SEPException(std::string("stream is not open correctly"));

    buf->resize(sz);

    stream.read(buf->getPtr(), sz);

    if (stream.received_hash() != stream.computed_hash())
      throw SEPException(std::string("Hashes do not match"));

    _bufferState = CPU_COMPRESSED;
  }
  if (_bufferState == UNDEFINED) throw SEPException("Bufferstate is undefined");
  return _buf->getSize() - oldSize;
}

gcpBuffer::gcpBuffer(
    const std::string &bucketName,
    google::cloud::v0::StatusOr<google::cloud::storage::Client> client,
    const std::string name, const std::vector<int> &n,
    const std::vector<int> &f, std::shared_ptr<compress> comp,
    const int ntrys) {
  setLoc(n, f);
  _client = client;
  _bucketName = bucketName;
  setName(name);
  setCompress(comp);
  setBufferState(ON_DISK);
  _ntrys = ntrys;
}
gcpBuffer::gcpBuffer(
    const std::string &bucketName,
    google::cloud::v0::StatusOr<google::cloud::storage::Client> client,
    const std::vector<int> &n, const std::vector<int> &f,
    std::shared_ptr<compress> comp, const bufferState state, const int ntrys) {
  _bucketName = bucketName;
  _client = client;

  setLoc(n, f);
  setCompress(comp);
  setBufferState(state);
  createStorageBuffer();
  _compress = comp;
  _ntrys = ntrys;
}
void gcpBuffer::remove() {
  namespace gcs = google::cloud::storage;

  _client.value().DeleteObject(_bucketName, _name);
}