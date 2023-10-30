
#include "document.h"
using namespace SEP;
document::document(std::shared_ptr<paramObj> par, std::string &commandLine,
                   std::string &shortDes) {
  _par = par;
  _com = commandLine;
  _short = shortDes;
}
void document::printDocument() {
  std::cerr << _com << "\n" << std::endl;
  std::cerr << "\t" << _short << std::endl;
  for (int i = 0; i < _pars.size(); i++) {
    std::cerr << _pars[i]._par << "\t" << _pars[i]._typ << std::endl;
    std::cerr << "\t [" << _pars[i]._val << "]\t" << _pars[i]._des;
  }
}

int document::getInt(const std::string &par, const std::string &des,
                     const int def) {
  _pars.push_back(paramO(par, std::string("int"), des));
  int x;
  try {
    x = _par->getInt(par, def);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(std::to_string(x));
  return x;
}
int document::getInt(const std::string &par, const std::string &des) {
  _pars.push_back(paramO(par, std::string("int"), des));
  int x;
  try {
    x = _par->getInt(par);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(std::to_string(x));

  return x;
}
float document::getFloat(const std::string &par, const std::string &des,
                         const float def) {
  _pars.push_back(paramO(par, std::string("float"), des));
  float x;
  try {
    x = _par->getFloat(par, def);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(std::to_string(x));

  return x;
}
float document::getFloat(const std::string &par, const std::string &des) {
  _pars.push_back(paramO(par, std::string("float"), des));
  float x;
  try {
    x = _par->getInt(par);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(std::to_string(x));

  return x;
}
std::string document::getString(const std::string &par, const std::string &des,
                                const std::string &def) {
  _pars.push_back(paramO(par, std::string("string"), des));

  std::string x;
  try {
    x = _par->getString(par, def);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(x);

  return x;
}
std::string document::getString(const std::string &par,
                                const std::string &des) {
  _pars.push_back(paramO(par, std::string("string"), des));
  std::string x;
  try {
    x = _par->getString(par);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(x);

  return x;
}
bool document::getBool(const std::string &par, const std::string &des,
                       const bool def) {
  _pars.push_back(paramO(par, std::string("bool"), des));

  bool x;
  try {
    x = _par->getBool(par, def);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(std::to_string(x));
  return x;
}
bool document::getBool(const std::string &par, const std::string &des) {
  _pars.push_back(paramO(par, std::string("bool"), des));
  bool x;
  try {
    x = _par->getBool(par);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  _pars[_pars.size() - 1].setValue(std::to_string(x));

  return x;
}
std::vector<int> document::getInts(const std::string &par,
                                   const std::string &des,
                                   const std::vector<int> def) {
  _pars.push_back(paramO(par, std::string("ints"), des));

  std::vector<int> x;
  try {
    x = _par->getInts(par, def);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  std::string m = std::to_string(x[0]);
  for (int i = 1; i < x.size(); i++) m = m + "," + std::to_string(x[i]);
  _pars[_pars.size() - 1].setValue(m);

  return x;
}
std::vector<int> document::getInts(const std::string &par, const int nval,
                                   const std::string &des) {
  _pars.push_back(paramO(par, std::string("ints"), des));

  std::vector<int> x;
  try {
    x = _par->getInts(par, nval);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  std::string m = std::to_string(x[0]);
  for (int i = 1; i < x.size(); i++) m = m + "," + std::to_string(x[i]);
  _pars[_pars.size() - 1].setValue(m);
  return x;
}
std::vector<float> document::getFloats(const std::string &par,
                                       const std::string &des,
                                       const std::vector<float> def) {
  _pars.push_back(paramO(par, std::string("floats"), des));

  std::vector<float> x;
  try {
    x = _par->getFloats(par, def);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  std::string m = std::to_string(x[0]);
  for (int i = 1; i < x.size(); i++) m = m + "," + std::to_string(x[i]);
  _pars[_pars.size() - 1].setValue(m);
  return x;
}
std::vector<float> document::getFloats(const std::string &par, const int nf,
                                       const std::string &des) {
  _pars.push_back(paramO(par, std::string("floats"), des));

  std::vector<float> x;
  try {
    x = _par->getFloats(par, nf);
  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  std::string m = std::to_string(x[0]);
  for (int i = 1; i < x.size(); i++) m = m + "," + std::to_string(x[i]);
  _pars[_pars.size() - 1].setValue(m);
  return x;
}

std::shared_ptr<genericRegFile> document::getRegFile(
    std::shared_ptr<genericIO> io, const std::string &par,
    const std::string &des, usage_code usage) {
  if (usage == usageIn)
    _pars.push_back(paramO(par, std::string("Input file"), des));
  else if (usage == usageOut)
    _pars.push_back(paramO(par, std::string("Output file"), des));
  else
    _pars.push_back(paramO(par, std::string("file"), des));

  std::shared_ptr<genericRegFile> x;
  try {
    x = io->getRegFile(par, usage);

  } catch (SEPException &x) {
    printDocument();
    std::cerr << x.what() << std::endl;
    exit(-1);
  }
  return x;
}
