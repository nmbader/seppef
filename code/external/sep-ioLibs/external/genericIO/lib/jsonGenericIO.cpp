#include "jsonGenericIO.h"
#include <exception>
#include <fstream>   // std::ifstream
#include <iostream>  // std::cout
using namespace SEP;
void jsonGenericIO::initJsonPars(const int argc, char **argv) {
  setValid(false);
  bool foundIn = false;
  std::string fileIn, fileOut;
  int i = 0;
  _progName = std::string(argv[0]);
  _type="JSON";
  /*Looking for jsonIn and jsonOut"*/
  while (i < argc && !foundIn) {
    std::string arg = std::string(argv[i]);
    if (arg.length() > 8) {
      if (arg.substr(0, 5) == std::string("json=")) {
        fileIn = arg.substr(5);
        foundIn = true;
      }
    }
    i++;
  }

  _init = false;
  if (!foundIn && !_sentError) {
    _sentError = true;

    _init = true;
    // inps.open(fileIn, std::ifstream::in);
    // std::shared_ptr<Json::Value> v(new Json::Value());
    // jsonArgs = v;

  } else {
    if (!inps) {
      std::cerr << std::string("Trouble opening " + std::string(fileIn))
                << std::endl;
      throw std::exception();
    }
    try {
      inps >> jsonArgs;
    } catch (int x) {
      std::cerr << std::string("Trouble parsing JSON file " +
                               std::string(fileIn))
                << std::endl;
      throw std::exception();
    }
  }
  std::shared_ptr<jsonParamObj> x(new jsonParamObj(jsonArgs));

  _type = "JSON";
  _param = x;
  setValid(true);
}
std::shared_ptr<genericRegFile> jsonGenericIO::getRegFileTag(
    const std::string &tag, const std::string &name, const usage_code usage,
    const int ndimMax) {
  if (!_init && !_sentError) {
    _sentError = true;
  }
  /*
     if((*jsonArgs)[name].isNull())  {
     std::cerr<<name<<std::string("  does not exist in json file")<<std::endl;
      throw std::exception();
     }
   */
  std::shared_ptr<jsonGenericRegFile> x(
      new jsonGenericRegFile(jsonArgs, usage, name, 0, 0, _progName, ndimMax));
  addRegFile(tag, x);
  return x;
}

std::shared_ptr<genericIrregFile> jsonGenericIO::getIrregFileTag(
    const std::string &tag, const std::string &name, const usage_code usage,
    const int ndimMax) {
  if (!_init && !_sentError) {
    _sentError = true;
  }

  std::shared_ptr<jsonGenericIrregFile> x(new jsonGenericIrregFile(
      jsonArgs, usage, name, 0, 0, _progName, ndimMax));
  addRegFile(tag, x);
  return x;
}

std::shared_ptr<paramObj> jsonGenericIO::getParamObj() {
  if (!_init && !_sentError) {
    _sentError = true;
  }
  return _param;
}
void jsonGenericIO::close() {
  for (auto i = _irregFiles.begin(); i != _irregFiles.end(); ++i) {
    std::shared_ptr<jsonGenericIrregFile> x =
        std::static_pointer_cast<jsonGenericIrregFile>(i->second);
    jsonArgs[i->first] = x->getArgs();
  }

  for (auto i = _regFiles.begin(); i != _regFiles.end(); ++i) {
    std::shared_ptr<jsonGenericRegFile> x =
        std::static_pointer_cast<jsonGenericRegFile>(i->second);
    jsonArgs[i->first] = x->getArgs();
  }
  filesClose();
}

/*
   std::static_pointer_cast<DerivedClass>(ptr_to_base)->f(); // OK
   // (constructs a temporary shared_ptr, then calls operator->)

   static_cast<DerivedClass*>(ptr_to_base.get())->f(); // also OK
   // (direct cast, does not construct a temporary shared_ptr)

 */
