#include "buffersIO.h"
#include <exception>
#include <fstream>   // std::ifstream
#include <iostream>  // std::cout
using namespace SEP;

std::shared_ptr<genericRegFile> buffersIO::getRegFileTag(
    const std::string &tag, const std::string &name, const usage_code usage) {
  if (!_init && !_sentError) {
    _sentError = true;
  }
  /*
     if((*jsonArgs)[name].isNull())  {
     std::cerr<<name<<std::string("  does not exist in json file")<<std::endl;
      throw std::exception();
     }
   */
  std::shared_ptr<jsonGenericFile> x(
      new buffersRegFile(jsonArgs, usage, name, _progName));
  addRegFile(tag, x);
  return x;
}
std::shared_ptr<genericIrregFile> buffersIO::getIrregFileTag(
    const std::string &tag, const std::string &name, const usage_code usage) {
  if (!_init && !_sentError) {
  }
  /*
     if((*jsonArgs)[name].isNull())  {
     std::cerr<<name<<std::string("  does not exist in json file")<<std::endl;
      throw std::exception();
     }
   */
  std::shared_ptr<jsonGenericFile> x(
      new buffersRegFile(jsonArgs, usage, name, _progName));
  addIrregFile(tag, x);
  return x;
}
std::shared_ptr<paramObj> buffersIO::getParamObj() {
  if (!_init && !_sentError) {
    std::cerr << std::string(
                     "json-file not specified no parameter reading possible")
              << std::endl;
    _sentError = true;
  }
  return _param;
}
void buffersIO::close() {
  for (auto i = _irregFiles.begin(); i != _irregFiles.end(); ++i) {
    std::shared_ptr<jsonGenericFile> x =
        std::static_pointer_cast<jsonGenericFile>(i->second);
    jsonArgs[i->first] = x->getArgs();
  }
  for (auto i = _regFiles.begin(); i != _regFiles.end(); ++i) {
    std::shared_ptr<jsonGenericFile> x =
        std::static_pointer_cast<jsonGenericFile>(i->second);
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
