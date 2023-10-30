#include "jsonGenericRegFile.h"
#include <cstdlib>
#include <exception>
#include<set>
#include <fstream>   // std::ifstream
#include <iostream>  // std::cout
using namespace SEP;
jsonGenericRegFile::jsonGenericRegFile(const Json::Value &arg,
                                       const usage_code usage,
                                       const std::string &tag, const int reelH,
                                       const int traceH,
                                       const std::string &progName,
                                       const int ndim) {
  _usage = usage;
  setupJson(arg, tag);
  _reelH = reelH;
  _traceH = traceH;

  if (!_newFile) {
    readDescription(ndim);
    std::shared_ptr<myFileIO> x(
        new myFileIO(getDataFileName(), usage, reelH, traceH,
                     jsonArgs.get("esize", 4).asInt(),
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = x;
  } else {
    std::string datapath = std::string("./");
    if (const char *env_p = std::getenv("DATAPATH"))
      datapath = std::string(env_p);
    _dataFile =
        datapath  + getJSONFileBaseName() + std::string(".dat");
    jsonArgs["filename"] = _dataFile;
  }
  jsonArgs["progName"] = progName;
}
void jsonGenericRegFile::setupJson(const Json::Value &arg,
                                   const std::string &tag,
                                   const std::string desFileDefault) {
  _tag = tag;

  if (arg[tag].isNull()) {
    _jsonFile = _tag + desFileDefault;
  } else {
    _jsonFile = arg[tag].asString();
  }
  _newFile = true;
  if (_usage == usageIn)
    _newFile = false;
  else if (_usage == usageInOut) {
    std::ifstream f(getJSONFileName());
    _newFile = !f.good();
    f.close();
  }

  if (_usage == usageIn || !_newFile) {
    std::ifstream inps;
    inps.open(getJSONFileName(), std::ifstream::in);
    if (!inps) {
      std::cerr << std::string("Trouble opening1 " + getJSONFileName())
                << std::endl;
      throw std::exception();
    }
    try {
      inps >> jsonArgs;
    } catch (int x) {
      std::cerr << std::string("Trouble parsing JSON file " + getJSONFileName())
                << std::endl;
      throw std::exception();
    }
    _dataFile = jsonArgs[std::string("filename")].asString();
  }
}

std::string jsonGenericRegFile::getJSONFileName() const { return _jsonFile; }


std::string jsonGenericRegFile::getJSONFileBaseName() const { 
  
  char seperator='/';
   std::string filePath=getJSONFileName();
	std::size_t sepPos = filePath.rfind(seperator);
  
	if(sepPos != std::string::npos)
	{
		return filePath.substr(sepPos + 1, filePath.size());
	}
	return filePath;
}
std::string jsonGenericRegFile::getDataFileName() const { return _dataFile; }
int jsonGenericRegFile::getInt(const std::string &arg) const {
  int x;
  if (jsonArgs[arg].isNull())
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  x = jsonArgs.get(arg, 1).asInt();
  return x;
}
int jsonGenericRegFile::getInt(const std::string &arg, const int def) const {
  int x = jsonArgs.get(arg, def).asInt();
  return x;
}
float jsonGenericRegFile::getFloat(const std::string &arg,
                                   const float def) const {
  float x;
  x = jsonArgs.get(arg, def).asFloat();
  return x;
}
float jsonGenericRegFile::getFloat(const std::string &arg) const {
  float x;
  if (jsonArgs[arg].isNull())
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  x = jsonArgs.get(arg, 1.).asFloat();
  return x;
}
std::string jsonGenericRegFile::getString(const std::string &arg) const {
  if (jsonArgs[arg].isNull())
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));

  return jsonArgs.get(arg, "").asString();
}
std::string jsonGenericRegFile::getString(const std::string &arg,
                                          const std::string &def) const {
  return jsonArgs.get(arg, def).asString();
}
bool jsonGenericRegFile::getBool(const std::string &arg, bool def) const {
  return jsonArgs.get(arg, def).asBool();
}
bool jsonGenericRegFile::getBool(const std::string &arg) const {
  if (jsonArgs[arg].isNull())
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  return jsonArgs.get(arg, false).asBool();
}
std::vector<int> jsonGenericRegFile::getInts(const std::string &arg,
                                             const int nvals) const {
  if (jsonArgs[arg].isNull())
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  const Json::Value vals = jsonArgs[arg];

  std::vector<int> x;
  for (int i = 0; i < nvals; i++) x.push_back(vals[i].asInt());
  return x;
}
std::vector<int> jsonGenericRegFile::getInts(
    const std::string &arg, const std::vector<int> &defs) const {
  std::vector<int> x;
  if (jsonArgs[arg].isNull()) {
    for (int i = 0; i < defs.size(); i++) x.push_back(defs[i]);
  } else {
    const Json::Value vals = jsonArgs[arg];
    for (int i = 0; i < defs.size(); i++) x.push_back(vals[i].asInt());
  }
  return x;
}
std::vector<float> jsonGenericRegFile::getFloats(const std::string &arg,
                                                 const int nvals) const {
  if (jsonArgs[arg].isNull())
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  const Json::Value vals = jsonArgs[arg];

  std::vector<float> x;
  for (int i = 0; i < nvals; i++) x.push_back(vals[i].asFloat());
  return x;
}
std::vector<float> jsonGenericRegFile::getFloats(
    const std::string &arg, const std::vector<float> &defs) const {
  std::vector<float> x;
  if (jsonArgs[arg].isNull()) {
    for (int i = 0; i < defs.size(); i++) x.push_back(defs[i]);
  } else {
    const Json::Value vals = jsonArgs[arg];
    for (int i = 0; i < defs.size(); i++) x.push_back(vals[i].asFloat());
  }
  return x;
}

void jsonGenericRegFile::putInt(const std::string &par, const int val) {
  jsonArgs[par] = val;
}
void jsonGenericRegFile::putFloat(const std::string &par, const float val) {
  jsonArgs[par] = val;
}
void jsonGenericRegFile::putString(const std::string &par,
                                   const std::string &val) {
  jsonArgs[par] = val;
}
void jsonGenericRegFile::putBool(const std::string &par, const bool val) {
  jsonArgs[par] = val;
}
void jsonGenericRegFile::putInts(const std::string &par,
                                 const std::vector<int> &val) {
  Json::Value vals;
  for (int i = 0; i < val.size(); i++) vals.append(val[i]);
  jsonArgs[par] = vals;
}
void jsonGenericRegFile::putFloats(const std::string &par,
                                   const std::vector<float> &val) {
  Json::Value vals;
  for (int i = 0; i < val.size(); i++) vals.append(val[i]);
  jsonArgs[par] = vals;
}
void jsonGenericRegFile::readDescription(const int ndimMax) {
  int ndim;
  bool breakIt = false;
  int iax = 9;
  while (iax >= 1 && !breakIt) {
    std::string tmp;
    int n = getInt("n" + std::to_string(iax), 1);
    float o = getFloat("o" + std::to_string(iax), 0.);

    if (n > 1 || fabsf(o) > 1e-4)
      breakIt = true;
    else
      iax--;
  }
  if (iax == 0) error("couldn't find any axes");
  ndim = iax;

  std::vector<axis> axes;
  for (int i = 1; i <= ndim; i++) {
    int n = getInt(std::string("n") + std::to_string(i), 1);
    float o = getFloat(std::string("o") + std::to_string(i), 0.);
    float d = getFloat(std::string("d") + std::to_string(i), 1.);
    std::string label = getString(std::string("label") + std::to_string(i), "");
    axes.push_back(axis(n, o, d, label));
  }
  std::string dtyp = getString(std::string("dataType"), std::string("FLOAT"));
  if (dtyp == std::string("FLOAT"))
    setDataType(DATA_FLOAT);
  else if (dtyp == std::string("COMPLEX"))
    setDataType(DATA_COMPLEX);
  else if (dtyp == std::string("INTEGER"))
    setDataType(DATA_INT);
  else if (dtyp == std::string("DOUBLE"))
    setDataType(DATA_DOUBLE);
      else if (dtyp == std::string("COMPLEXDOUBLE"))
    setDataType(DATA_COMPLEXDOUBLE);
  else if (dtyp == std::string("BYTE"))
    setDataType(DATA_BYTE);

  std::shared_ptr<hypercube> hyper(new hypercube(axes));
  setHyper(hyper);
}
void jsonGenericRegFile::setHistory(const Json::Value &input) {
  Json::Value hist;
  if (!input["history"].isNull()) {
    hist = input["history"];
  }
  Json::Value last;
  for (Json::Value::const_iterator itr = input.begin(); itr != input.end();
       itr++) {
    std::string key = itr.key().toStyledString();
    Json::Value value = (*itr);
    if (key == std::string("history")) {
      last[key] = value;
    }
  }
  hist.append(last);
  jsonArgs["history"] = hist;
}

void jsonGenericRegFile::writeDescription() {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<axis> axes = hyper->returnAxes(hyper->getNdim());
  for (int i = 1; i <= axes.size(); i++) {
    putInt(std::string("n") + std::to_string(i), axes[i - 1].n);
    putFloat(std::string("o") + std::to_string(i), axes[i - 1].o);
    putFloat(std::string("d") + std::to_string(i), axes[i - 1].d);
    putString(std::string("label") + std::to_string(i), axes[i - 1].label);
  }
}
void jsonGenericRegFile::close() {
  myio->close();
  if (_usage == usageOut || _usage == usageInOut) {
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
  }
}
void jsonGenericRegFile::readFloatStream(float *array, const long long npts) {
  long long maxsize = 10000000;
  long long nread = 0;
  long long nptsT = npts * 4;
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }

  myio->readTraceStream(npts, array);
}
void jsonGenericRegFile::readDoubleStream(double *array, const long long npts) {
  long long maxsize = 10000000;
  long long nread = 0;
  long long nptsT = npts * 8;
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }

  myio->readTraceStream(npts, array);
}

void jsonGenericRegFile::readIntStream(int *array, const long long npts) {
  long long maxsize = 10000000;
  long long nread = 0;
  long long nptsT = npts * 4;
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }

  myio->readTraceStream(npts, array);
}

void jsonGenericRegFile::readByteStream(unsigned char *array,
                                        const long long npts) {
  long long maxsize = 10000000;
  long long nread = 0;
  long long nptsT = npts * 8;
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }

  myio->readTraceStream(npts, array);
}

void jsonGenericRegFile::seekTo(const long long iv, const int whence) {
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->seekTo(iv, whence);
}

void jsonGenericRegFile::writeFloatStream(const float *array,
                                          const long long npts) {
  long long maxsize = 10000000;
  long long nwrite = 0;
  long long nptsT = npts * 4;
  setDataType(DATA_FLOAT);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeTraceStream(npts, array);
}

void jsonGenericRegFile::writeByteStream(const unsigned char *array,
                                         const long long npts) {
  long long maxsize = 10000000;
  long long nwrite = 0;
  long long nptsT = npts * 1;
  setDataType(DATA_FLOAT);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 1,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeTraceStream(npts, array);
}

void jsonGenericRegFile::writeDoubleStream(const double *array,
                                           const long long npts) {
  long long maxsize = 10000000;
  long long nwrite = 0;
  long long nptsT = npts * 8;
  setDataType(DATA_DOUBLE);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeTraceStream(npts, array);
}

void jsonGenericRegFile::writeIntStream(const int *array,
                                        const long long npts) {
  long long maxsize = 10000000;
  long long nwrite = 0;
  long long nptsT = npts * 4;
  setDataType(DATA_INT);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeTraceStream(npts, array);
}

void jsonGenericRegFile::readFloatWindow(const std::vector<int> &nw,
                                         const std::vector<int> &fw,
                                         const std::vector<int> &jw,
                                         float *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->readWindow(nw, fw, jw, array);
}

void jsonGenericRegFile::readDoubleWindow(const std::vector<int> &nw,
                                          const std::vector<int> &fw,
                                          const std::vector<int> &jw,
                                          double *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->readWindow(nw, fw, jw, array);
}

void jsonGenericRegFile::readIntWindow(const std::vector<int> &nw,
                                       const std::vector<int> &fw,
                                       const std::vector<int> &jw, int *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->readWindow(nw, fw, jw, array);
}

long long jsonGenericRegFile::getDataSize() {
  if (!myio) {
    std::shared_ptr<myFileIO> iox(new myFileIO(
        getDataFileName(), _usage, _reelH, _traceH, 1, false, getHyper()));
    myio = iox;
  }
  return myio->getSize();
}
void jsonGenericRegFile::readByteWindow(const std::vector<int> &nw,
                                        const std::vector<int> &fw,
                                        const std::vector<int> &jw,
                                        unsigned char *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(new myFileIO(
        getDataFileName(), _usage, _reelH, _traceH, 1, false, getHyper()));
    myio = iox;
  }
  myio->readWindow(nw, fw, jw, array);
}
void jsonGenericRegFile::writeFloatWindow(const std::vector<int> &nw,
                                          const std::vector<int> &fw,
                                          const std::vector<int> &jw,
                                          const float *array) {
  setDataType(DATA_FLOAT);




  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();

  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }

  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }

  if (!myio) {

    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    
    myio = iox;
  }

  myio->writeWindow(nw, fw, jw, array);
}

void jsonGenericRegFile::writeByteWindow(const std::vector<int> &nw,
                                         const std::vector<int> &fw,
                                         const std::vector<int> &jw,
                                         const unsigned char *array) {
  setDataType(DATA_BYTE);

  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 1,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeWindow(nw, fw, jw, array);
}

void jsonGenericRegFile::writeIntWindow(const std::vector<int> &nw,
                                        const std::vector<int> &fw,
                                        const std::vector<int> &jw,
                                        const int *array) {
  setDataType(DATA_INT);

  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeWindow(nw, fw, jw, array);
}

void jsonGenericRegFile::writeDoubleWindow(const std::vector<int> &nw,
                                           const std::vector<int> &fw,
                                           const std::vector<int> &jw,
                                           const double *array) {
  setDataType(DATA_DOUBLE);

  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 4,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeWindow(nw, fw, jw, array);
}

void jsonGenericRegFile::readComplexStream(std::complex<float> *array,
                                           const long long npts) {
  long long maxsize = 10000000;
  long long nread = 0;
  long long nptsT = npts * 4;
  setDataType(DATA_COMPLEX);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->readTraceStream(npts, array);
}

void jsonGenericRegFile::readComplexDoubleStream(std::complex<double> *array,
                                           const long long npts) {
  long long maxsize = 10000000;
  long long nread = 0;
  long long nptsT = npts * 4;
  setDataType(DATA_COMPLEXDOUBLE);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 16,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->readTraceStream(npts, array);
}
void jsonGenericRegFile::writeComplexStream(const std::complex<float> *array,
                                            const long long npts) {
  long long maxsize = 10000000;
  long long nwrite = 0;
  long long nptsT = npts * 4;
  setDataType(DATA_COMPLEX);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeTraceStream(npts, array);
}
void jsonGenericRegFile::writeComplexDoubleStream(const std::complex<double> *array,
                                            const long long npts) {
  long long maxsize = 10000000;
  long long nwrite = 0;
  long long nptsT = npts * 4;
  setDataType(DATA_COMPLEXDOUBLE);

  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 16,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeTraceStream(npts, array);
}

void jsonGenericRegFile::readComplexWindow(const std::vector<int> &nw,
                                           const std::vector<int> &fw,
                                           const std::vector<int> &jw,
                                           std::complex<float> *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  setDataType(DATA_COMPLEX);

  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->readWindow(nw, fw, jw, array);
}

void jsonGenericRegFile::readComplexDoubleWindow(const std::vector<int> &nw,
                                           const std::vector<int> &fw,
                                           const std::vector<int> &jw,
                                           std::complex<double> *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  setDataType(DATA_COMPLEXDOUBLE);

  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 16,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->readWindow(nw, fw, jw, array);
}
void jsonGenericRegFile::writeComplexWindow(const std::vector<int> &nw,
                                            const std::vector<int> &fw,
                                            const std::vector<int> &jw,
                                            const std::complex<float> *array) {
  setDataType(DATA_COMPLEX);

  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 8,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeWindow(nw, fw, jw, array);
}
void jsonGenericRegFile::writeComplexDoubleWindow(const std::vector<int> &nw,
                                            const std::vector<int> &fw,
                                            const std::vector<int> &jw,
                                            const std::complex<double> *array) {
  setDataType(DATA_COMPLEX);

  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }
  if (!myio) {
    std::shared_ptr<myFileIO> iox(
        new myFileIO(getDataFileName(), _usage, _reelH, _traceH, 16,
                     jsonArgs.get("swapData", false).asBool(), getHyper()));
    myio = iox;
  }
  myio->writeWindow(nw, fw, jw, array);
}
void jsonGenericRegFile::message(const std::string &errm) const {
  std::cerr << errm << std::endl;
}
void jsonGenericRegFile::error(const std::string &errm) const {
  std::cerr << errm << std::endl;
  throw std::exception();
}
