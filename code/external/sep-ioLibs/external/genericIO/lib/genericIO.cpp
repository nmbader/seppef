#include "genericIO.h"
using namespace SEP;
void genericIO::filesClose() {
  for (auto i = _regFiles.begin(); i != _regFiles.end(); ++i) {
    i->second->close();
  }
  for (auto i = _irregFiles.begin(); i != _irregFiles.end(); ++i) {
    i->second->close();
  }
}
std::shared_ptr<SEP::genericRegFile> genericIO::getRegFile(
    const std::string& name, const SEP::usage_code usage, const int ndimMax) {
  std::shared_ptr<paramObj> par = getParamObj();
  std::string filename = par->getString(name, name);
  return getRegFileTag(name, filename, usage, ndimMax);
}

std::shared_ptr<SEP::genericRegFile> genericIO::getDocRegFile(
    const std::string& name, const std::string& doc,
    const SEP::usage_code usage, const int ndimMax) {
  std::shared_ptr<paramObj> par = getParamObj();
  std::string filename = par->getString(name, name);
  std::shared_ptr<SEP::genericRegFile> v;
  try {
    v = getRegFileTag(name, filename, usage, ndimMax);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}

void genericIO::removeRegFile(const std::string fle) {
  if (_regFiles.count(fle) == 0)
    throw SEPException(std::string("Requested to remove file ") + fle +
                       std::string("that has not been initialized"));
  _regFiles[fle]->remove();
  _regFiles.erase(fle);
}

std::shared_ptr<SEP::genericRegFile> genericIO::getDocRegFile(
    const std::string& name, const std::string& doc, const std::string usage,
    const int ndimMax) {
  std::shared_ptr<SEP::genericRegFile> v;
  try {
    v = getRegFile(name, usage, ndimMax);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}

std::shared_ptr<SEP::genericIrregFile> genericIO::getIrregFile(
    const std::string& name, const SEP::usage_code usage, const int ndimMax) {
  std::shared_ptr<paramObj> par = getParamObj();
  std::string filename = par->getString(name, name);

  return getIrregFileTag(name, filename, usage, ndimMax);
}
void genericIO::fileDebug(const std::string nm, const float* data, const int n1,
                          const int n2, const int n3) {
  std::vector<SEP::axis> axs;
  axs.push_back(SEP::axis(n1));
  axs.push_back(SEP::axis(n2));
  axs.push_back(SEP::axis(n3));
  std::shared_ptr<hypercube> hyper2(new hypercube(axs));

  std::shared_ptr<SEP::genericRegFile> fle = getRegFile(nm, SEP::usageOut);
  fle->setHyper(hyper2);
  fle->writeDescription();
  fle->writeFloatStream(data, (long long)n1 * (long long)n2 * (long long)n3);
}
void genericIO::fileDebug(const std::string nm, const float* data, const int n1,
                          const int n2) {
  std::vector<SEP::axis> axs;
  axs.push_back(SEP::axis(n1));
  axs.push_back(SEP::axis(n2));
  std::shared_ptr<hypercube> hyper2(new hypercube(axs));

  std::shared_ptr<SEP::genericRegFile> fle = getRegFile(nm, SEP::usageOut);
  fle->setHyper(hyper2);
  fle->writeDescription();
  fle->writeFloatStream(data, (long long)n1 * (long long)n2);
}

void genericIO::addFileDescription(const std::shared_ptr<genericRegFile> fileIn,
                                   const std::string& title,
                                   std::shared_ptr<genericRegFile> fileOut) {
  Json::Value val = fileIn->getDescription();
  fileOut->putDescription(title, val);
}
