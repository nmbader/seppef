#include "paramObj.h"
using namespace SEP;

int paramObj::getDocInt(const std::string& arg, const std::string& doc) {
  int v;
  try {
    v = getInt(arg);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}

int paramObj::getDocInt(const std::string& arg, const std::string& doc,
                        const int def) {
  int v;
  try {
    v = getInt(arg, def);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
float paramObj::getDocFloat(const std::string& arg, const std::string& doc,
                            const float def) {
  float v;
  try {
    v = getFloat(arg, def);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
float paramObj::getDocFloat(const std::string& arg, const std::string& doc) {
  float v;
  try {
    v = getFloat(arg);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
std::string paramObj::getDocString(const std::string& arg,
                                   const std::string& doc) {
  std::string v;
  try {
    v = getString(arg);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
std::string paramObj::getDocString(const std::string& arg,
                                   const std::string& doc,
                                   const std::string& def) {
  std::string v;
  try {
    v = getString(arg, def);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
bool paramObj::getDocBool(const std::string& arg, const std::string& doc,
                          const bool def) {
  bool v;
  try {
    v = getBool(arg, def);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
bool paramObj::getDocBool(const std::string& arg, const std::string& doc) {
  bool v;
  try {
    v = getBool(arg);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
std::vector<int> paramObj::getDocInts(const std::string& arg,
                                      const std::string& doc, const int nvals) {
  std::vector<int> v;
  try {
    v = getInts(arg, nvals);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
std::vector<int> paramObj::getDocInts(const std::string& arg,
                                      const std::string& doc,
                                      const std::vector<int>& defs) {
  std::vector<int> v;
  try {
    v = getInts(arg, defs);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
std::vector<float> paramObj::getDocFloats(const std::string& arg,
                                          const std::string& doc, int nvals) {
  std::vector<float> v;
  try {
    v = getFloats(arg, nvals);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
std::vector<float> paramObj::getDocFloats(const std::string& arg,
                                          const std::string& doc,
                                          const std::vector<float>& defs) {
  std::vector<float> v;
  try {
    v = getFloats(arg, defs);
  } catch (SEPException& x) {
    std::string tmp = std::string(x.what()) + std::string("\n\t") + doc;
    throw SEPException(tmp);
  }
  return v;
}
void paramObj::addParams(std::map<std::string, std::string>& pars) {
  throw SEPException(std::string("Must override addParams"));
}
void paramObj::addParams(std::vector<std::string>& pars) {
  for (auto a : pars) {
    int index = a.find("=");
    if (index != std::string::npos) {
      std::map<std::string, std::string> x;
      x[a.substr(0, index)] = a.substr(index + 1);
      addParams(x);
    }
  }
}
