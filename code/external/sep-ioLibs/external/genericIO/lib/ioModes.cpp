#include "ioModes.h"
#include <iostream>
#include "buffersConfig.h"
#include "dictParams.h"
#include "fileBuffersIO.h"
#include "ioConfig.h"
#include "memoryIO.h"
#include "segyIO.h"
#ifdef USE_RSF
#include "rsfIO.h"
#endif
#ifdef USE_SEP
#include "sepIO.h"
#endif
#ifdef USE_GCP
#include "gcpBuffersIO.h"
#endif
using namespace SEP;
void ioModes::setup(const int argc, char **argv) {
  std::shared_ptr<jsonGenericIO> a(new jsonGenericIO(argc, argv));

  _ios["JSON"] = a;

  std::shared_ptr<segyIO> d(new segyIO(argc, argv));
  _ios["SEGY"] = d;

  std::shared_ptr<fileBuffersIO> m(new fileBuffersIO(argc, argv));
  _ios["FILEBUFFERS"] = m;

#ifdef USE_GCP
  std::shared_ptr<gcpBuffersIO> y(new gcpBuffersIO(argc, argv));
  _ios["GCPBUFFERS"] = y;
#endif

#ifdef USE_RSF
  std::shared_ptr<rsfIO> b(new rsfIO(argc, argv));
  _ios["RSF"] = b;
#endif
#ifdef USE_SEP

  std::shared_ptr<sepIO> c(new sepIO(argc, argv));
  _ios["SEP"] = c;
#else
#endif
  std::map<std::string, std::string> dict;
  std::shared_ptr<memoryIO> e(new memoryIO(dict));
  _ios["memory"] = e;

  _defaultType = DEFAULTIO;
  _defaultIO = _ios[_defaultType];

  std::shared_ptr<paramObj> _par = getParamObj();
  changeParamObj(_par);
}
std::shared_ptr<genericIO> ioModes::getDefaultIO() {
  return getIO(_defaultType);
}

std::shared_ptr<genericIO> ioModes::getIO(const std::string &def) {
  if (_ios.count(def) != 1)
    throw SEPException(def +
                       std::string(" io has not been defined and/or built"));
  return _ios[def];
}
std::shared_ptr<genericRegFile> ioModes::getRegFileTag(
    const std::string &tag, const std::string &ioname, const std::string &name,
    usage_code usage) {
  if (_ios.count(ioname) == 0)
    _par->error(ioname + " io has not been ioname and/or built");
  if (!_ios[ioname]->getValid())
    _par->error(ioname + std::string(" has not been initialized correctly"));
  return _ios[ioname]->getRegFile(name, usage);
}
std::shared_ptr<genericRegFile> ioModes::getGenericRegFile(
    const std::string &name, const usage_code usage) {
  return _defaultIO->getRegFile(name, usage);
}

std::shared_ptr<genericIO> ioModes::getInputIO() {
  std::string iouse = getParamObj()->getString(
      "inputIO", getParamObj()->getString("IO", std::string("DEFAULT")));
  if (strcmp(iouse.c_str(), "DEFAULT") == 0) return getDefaultIO();
  try {
    return getIO(iouse);
  } catch (SEPException &x) {
    throw x;
  }
}
std::vector<std::string> ioModes::getIOs() const {
  std::vector<std::string> ios;
  for (auto io = _ios.begin(); io != _ios.end(); ++io) {
    ios.push_back(io->first);
  }
  return ios;
}
std::shared_ptr<genericIO> ioModes::getOutputIO() {
  std::string iouse = getParamObj()->getString(
      "outputIO", getParamObj()->getString("IO", std::string("DEFAULT")));
  if (strcmp(iouse.c_str(), "DEFAULT") == 0) return getDefaultIO();
  try {
    return getIO(iouse);
  } catch (SEPException &x) {
    throw x;
  }
}
void ioModes::changeParamObj(std::shared_ptr<paramObj> par) {
  for (auto it = _ios.begin(); it != _ios.end(); it++) {
    it->second->setParamObj(par);
  }
}
std::shared_ptr<paramObj> ioModes::getParamObj() {
  std::string iouse = getDefaultIO()->getParamObj()->getString(
      "paramIO",
      getDefaultIO()->getParamObj()->getString("IO", std::string("DEFAULT")));
  if (strcmp(iouse.c_str(), "DEFAULT") == 0)
    return getDefaultIO()->getParamObj();
  try {
    return getIO(iouse)->getParamObj();
  } catch (SEPException &x) {
    throw x;
  }
}

std::shared_ptr<ioModesFortran> ioModesFortran::instance = nullptr;
void ioModesFortran::setup(const int argc, char **argv) {
  std::shared_ptr<ioModes> x(new ioModes(argc, argv));
  _io = x;
}
