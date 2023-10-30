#include "fileBuffersRegFile.h"
#include <unistd.h>
#include <cstdlib>
#include <exception>
#include <fstream>   // std::ifstream
#include <iostream>  // std::cout
#include "fileBuffers.h"
using namespace SEP;
fileBuffersRegFile::fileBuffersRegFile(const Json::Value &arg,
                                       const usage_code usage,
                                       const std::string &tag,
                                       const std::string &progName,
                                       const int ndimMax) {
  setUsage(usage);
  setupJson(arg, tag, std::string("/des.dat"));

  if (!_newFile) {
    readDescription(ndimMax);

    if (jsonArgs["bufferInfo"].isNull())
      error(std::string("bufferInfo not provided in JSON file"));
    _bufs.reset(
        new SEP::IO::fileBuffers(getHyper(), tag, jsonArgs["bufferInfo"]));
  }
  _usage = usage;
  jsonArgs["progName"] = progName;
  jsonArgs["name"] = tag;
}

void fileBuffersRegFile::close() {
  if (getUsage() == usageOut || getUsage() == usageInOut) {
    std::ofstream outps;
    outps.open(getJSONFileName(), std::ofstream::out);
    if (!outps) {
      std::cerr << std::string("Trouble opening for write") + getJSONFileName()
                << std::endl;
      throw std::exception();
    }
    try {
      outps << jsonArgs;
    } catch (int x) {
      std::cerr << std::string("Trouble writing JSON file ") + getJSONFileName()
                << std::endl;
      throw std::exception();
    }
    _bufs->changeState(SEP::IO::ON_DISK);
  }
}
void fileBuffersRegFile::createBuffers() {
  if (_bufs) return;
  if (!_hyper) error("Must set hypercube before blocking");
  if (getDataType() == SEP::DATA_UNKNOWN)
    error("Must set dataType before setting blocks");
  _bufs.reset(
      new SEP::IO::fileBuffers(getHyper(), getDataType(), _comp, _block, _mem));
  _bufs->setName(jsonArgs["name"].asString(), true);
}

void fileBuffersRegFile::removeDescDir() {
  std::string desc = _tag + std::string("/des.dat");
  std::remove(desc.c_str());
  rmdir(_tag.c_str());
}

fileBuffersIrregFile::fileBuffersIrregFile(const Json::Value &arg,
                                           const usage_code usage,
                                           const std::string &tag,
                                           const std::string &progName,
                                           const int ndimMax) {
  setUsage(usage);
  setupJson(arg, tag, std::string("/des.dat"));

  if (!_newFile) {
    readDescription(ndimMax);

    if (jsonArgs["bufferInfo"].isNull())
      error(std::string("bufferInfo not provided in JSON file"));
    _bufs.reset(
        new SEP::IO::fileBuffers(getHyper(), tag, jsonArgs["bufferInfo"]));
  }

  jsonArgs["progName"] = progName;
  jsonArgs["name"] = tag;
}

void fileBuffersIrregFile::close() {
  if (getUsage() == usageOut || getUsage() == usageInOut) {
    std::ofstream outps;
    outps.open(getJSONFileName(), std::ofstream::out);
    if (!outps) {
      std::cerr << std::string("Trouble opening for write") + getJSONFileName()
                << std::endl;
      throw std::exception();
    }
    try {
      outps << jsonArgs;
    } catch (int x) {
      std::cerr << std::string("Trouble writing JSON file ") + getJSONFileName()
                << std::endl;
      throw std::exception();
    }
    _bufs->changeState(SEP::IO::ON_DISK);
  }
}
void fileBuffersIrregFile::createBuffers() {
  if (_bufs) return;
  if (!_hyper) error("Must set hypercube before blocking");
  if (getDataType() == SEP::DATA_UNKNOWN)
    error("Must set dataType before setting blocks");
  _bufs.reset(
      new SEP::IO::fileBuffers(getHyper(), getDataType(), _comp, _block, _mem));
  _bufs->setName(jsonArgs["name"].asString(), true);
}
