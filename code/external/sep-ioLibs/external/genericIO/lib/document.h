#ifndef DOCUMENT_H
#define DOCUMENT_H 1
#include <genericIO.h>
#include <paramObj.h>
#include <string>

namespace SEP {

class paramO {
 public:
  paramO(const std::string &par, const std::string &typ,
         const std::string &des) {
    _par = par;
    _typ = typ;
    _des = des;
    _val = std::string("");
  }
  void setValue(const std::string val) { _val = val; }

  std::string _par;  //< Parmeter name
  std::string _typ;  //< Parameter type
  std::string _des;  //< Parameter description
  std::string _val;  //< Parameter value
};

class document {
 public:
  document(std::shared_ptr<paramObj> par, std::string &commandLine,
           std::string &shortDes);
  void printDocument();
  int getInt(const std::string &par, const std::string &des, const int def);
  int getInt(const std::string &par, const std::string &des);
  float getFloat(const std::string &par, const std::string &des,
                 const float def);
  float getFloat(const std::string &par, const std::string &des);
  std::string getString(const std::string &par, const std::string &des,
                        const std::string &def);
  std::string getString(const std::string &par, const std::string &des);
  bool getBool(const std::string &par, const std::string &des);
  bool getBool(const std::string &par, const std::string &des, const bool def);
  std::vector<int> getInts(const std::string &par, const std::string &des,
                           const std::vector<int> def);
  std::vector<int> getInts(const std::string &par, const int nf,
                           const std::string &des);
  std::vector<float> getFloats(const std::string &par, const std::string &des,
                               const std::vector<float> def);
  std::vector<float> getFloats(const std::string &par, const int nf,
                               const std::string &des);

  std::shared_ptr<genericRegFile> getRegFile(std::shared_ptr<genericIO> io,
                                             const std::string &par,
                                             const std::string &des,
                                             SEP::usage_code usage);

 private:
  std::shared_ptr<paramObj> _par;
  std::vector<paramO> _pars;
  std::string _com, _short;
};

}  // namespace SEP
#endif