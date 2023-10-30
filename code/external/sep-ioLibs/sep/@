#include "sep_reg_file.h"
extern "C" {
#include "sep3d.h"
#include "seplib.h"
}
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace SEP;
sepRegFile::sepRegFile(const std::string &tag, const usage_code usage,
                       const int ndim) {
  _tag = tag;
  _usage = usage;
  switch (usage) {
    case usageIn:

      if (_tag != "in")
        if (NULL == auxin(_tag.c_str()))
          error(std::string("can not open file ") + tag);

      readDescription(ndim);

      break;
    case usageOut:
      if (tag != "out")
        if (0 == auxout(_tag.c_str()))
          error(std::string("can not open file ") + tag);
      break;
    case usageInOut:
      if (0 == auxinout(_tag.c_str()))
        error(std::string("can not open file ") + tag);
      break;
    case usageScr:
      if (0 == auxscr(_tag.c_str()))
        error(std::string("can not open file ") + tag);
      break;
    default:
      error("can't handle type");
  }
}

void sepRegFile::remove() {
  char temp[20000];
  auxpar("in", "s", temp, _tag.c_str());
  std::stringstream test(temp);
  std::string segment;

  while (std::getline(test, segment, ';'))
    std::remove(segment.c_str());  // delete file

  strcpy(temp, _tag.c_str());
  auxpar(_tag.c_str(), "s", temp, _tag.c_str());
  segment = temp;

  std::remove(segment.c_str());  // delete file
}
int sepRegFile::getInt(const std::string &arg) const {
  int x;
  if (0 == auxpar(arg.c_str(), "d", &x, _tag.c_str()))
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  return x;
}
int sepRegFile::getInt(const std::string &arg, const int def) const {
  int x = def;
  int i = auxpar(arg.c_str(), "d", &x, _tag.c_str());
  return x;
}

float sepRegFile::getFloat(const std::string &arg, const float def) const {
  float x = def;
  int i = auxpar(arg.c_str(), "f", &x, _tag.c_str());
  return x;
}
float sepRegFile::getFloat(const std::string &arg) const {
  float x;
  if (0 == auxpar(arg.c_str(), "f", &x, _tag.c_str()))
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  return x;
}

std::string sepRegFile::getString(const std::string &arg) const {
  char buf[10000];
  if (0 == auxpar(arg.c_str(), "s", buf, _tag.c_str()))
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  return std::string(buf);
}
std::string sepRegFile::getString(const std::string &arg,
                                  const std::string &def) const {
  char buf[10000];
  std::copy(def.begin(), def.end(), buf);
  int i = auxpar(arg.c_str(), "s", buf, _tag.c_str());
  return std::string(buf);
}

bool sepRegFile::getBool(const std::string &arg, const bool def) const {
  bool x = def;
  int i = auxpar(arg.c_str(), "l", &x, _tag.c_str());
  return x;
}
bool sepRegFile::getBool(const std::string &arg) const {
  bool x;
  if (0 == auxpar(arg.c_str(), "l", &x, _tag.c_str())) {
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  }
  return x;
}

std::vector<int> sepRegFile::getInts(const std::string &arg, int num) const {
  int tmp[10000];
  int ierr = auxpar(arg.c_str(), "d", tmp, _tag.c_str());
  if (ierr == 0)
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  std::vector<int> x;
  for (int i = 0; i < ierr; i++) x.push_back(tmp[i]);
  return x;
}
std::vector<int> sepRegFile::getInts(const std::string &arg,
                                     const std::vector<int> &defs) const {
  int tmp[10000];
  for (int i = 0; i < defs.size(); i++) {
    tmp[i] = defs[i];
  }
  int ierr = auxpar(arg.c_str(), "d", tmp, _tag.c_str());
  if (ierr == 0)
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  std::vector<int> x;
  if (ierr > 0) {
    for (int i = 0; i < ierr; i++) x.push_back(tmp[i]);
  } else {
    for (int i = 0; i < defs.size(); i++) x.push_back(defs[i]);
  }
  return x;
}

std::vector<float> sepRegFile::getFloats(const std::string &arg,
                                         int num) const {
  float tmp[10000];
  int ierr = auxpar(arg.c_str(), "f", tmp, _tag.c_str());
  if (ierr == 0)
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  std::vector<float> x;
  for (int i = 0; i < ierr; i++) x.push_back(tmp[i]);
  return x;
}
std::vector<float> sepRegFile::getFloats(const std::string &arg,
                                         const std::vector<float> &defs) const {
  float tmp[10000];
  for (int i = 0; i < defs.size(); i++) {
    tmp[i] = defs[i];
  }
  int ierr = auxpar(arg.c_str(), "f", tmp, _tag.c_str());
  if (ierr == 0)
    error(std::string("trouble grabbing parameter ") + arg +
          std::string(" from parameters"));
  std::vector<float> x;
  if (ierr > 0) {
    for (int i = 0; i < ierr; i++) x.push_back(tmp[i]);
  } else {
    for (int i = 0; i < defs.size(); i++) x.push_back(defs[i]);
  }
  return x;
}
void sepRegFile::close() {
  if (_tag == std::string("in")) {
    hclose();
  } else
    auxclose(_tag.c_str());
}
void sepRegFile::error(const std::string &err) const { seperr(err.c_str()); }

void sepRegFile::putInt(const std::string &par, const int val) {
  auxputch(par.c_str(), "d", &val, _tag.c_str());
}
void sepRegFile::putFloat(const std::string &par, const float val) {
  auxputch(par.c_str(), "f", &val, _tag.c_str());
}
void sepRegFile::putString(const std::string &par, const std::string &val) {
  char x[9999];
  std::copy(val.begin(), val.end(), x);
  x[val.length()] = '\0';
  auxputch(par.c_str(), "s", val.c_str(), _tag.c_str());
}

void sepRegFile::putBool(const std::string &par, const bool val) {
  int x = 0;
  if (val) x = 1;
  auxputch(par.c_str(), "l", &x, _tag.c_str());
}
void sepRegFile::putInts(const std::string &par, const std::vector<int> &val) {
  int *tmp = new int[val.size()];
  for (int i = 0; i < val.size(); i++) tmp[i] = val[i];
  auxputch(par.c_str(), "d", tmp, _tag.c_str());
  delete[] tmp;
}
void sepRegFile::putFloats(const std::string &par,
                           const std::vector<float> &val) {
  float *tmp = new float[val.size()];
  for (int i = 0; i < val.size(); i++) tmp[i] = val[i];
  auxputch(par.c_str(), "f", tmp, _tag.c_str());
  delete[] tmp;
}
void sepRegFile::readFloatStream(float *array, const long long npts) {
  long long nptsT = npts * 4;
  long long ierr = sreed_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble reading from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}
void sepRegFile::readIntStream(int *array, const long long npts) {
  long long nptsT = npts * 4;
  long long ierr = sreed_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble reading from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}
void sepRegFile::readDoubleStream(double

                                      *array,
                                  const long long npts) {
  long long nptsT = npts * 8;
  long long ierr = sreed_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble reading from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}
void sepRegFile::readComplexStream(std::complex<float> *array,
                                   const long long npts) {
  long long nptsT = npts * 8;
  long long ierr = sreed_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble reading from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}
void sepRegFile::readComplexDoubleStream(std::complex<double> *array,
                                   const long long npts) {
  long long nptsT = npts * 16;
  long long ierr = sreed_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble reading from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}

void sepRegFile::readByteStream(unsigned char *array, const long long npts) {
  long long nptsT = npts * 1;
  long long ierr = sreed_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble reading from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}

void sepRegFile::writeFloatStream(const float *array, const long long npts) {
  long long nptsT = npts * 4;
  long long ierr = srite_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble write from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}
void sepRegFile::writeByteStream(const unsigned char *array,
                                 const long long npts) {
  long long nptsT = npts * 1;
  long long ierr = srite_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble write from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}
void sepRegFile::writeIntStream(const int *array, const long long npts) {
  long long nptsT = npts * 4;
  long long ierr = srite_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble write from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}
void sepRegFile::writeDoubleStream(const double *array, const long long npts) {
  long long nptsT = npts * 8;
  set_format(_tag.c_str(), "native_double");
  long long ierr = srite_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble write from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}

void sepRegFile::readFloatWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw, float *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();

  setDataType(DATA_FLOAT);

  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }

  int ndim = ng.size();
  if (0 != sreed_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 4, array))
    error(std::string("trouble reading data from tag ") + _tag);
}

void sepRegFile::readDoubleWindow(const std::vector<int> &nw,
                                  const std::vector<int> &fw,
                                  const std::vector<int> &jw, double *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();

  setDataType(DATA_DOUBLE);

  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }

  int ndim = ng.size();
  if (0 != sreed_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 8, array))
    error(std::string("trouble reading data from tag ") + _tag);
}
void sepRegFile::readIntWindow(const std::vector<int> &nw,
                               const std::vector<int> &fw,
                               const std::vector<int> &jw, int *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<int> ng = hyper->getNs();

  setDataType(DATA_INT);

  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }

  int ndim = ng.size();
  if (0 != sreed_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 4, array))
    error(std::string("trouble reading data from tag ") + _tag);
}

void sepRegFile::readComplexWindow(const std::vector<int> &nw,
                                   const std::vector<int> &fw,
                                   const std::vector<int> &jw,
                                   std::complex<float> *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  setDataType(DATA_COMPLEX);

  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }

  int ndim = ng.size();
  if (0 != sreed_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 8, array))
    error(std::string("trouble reading data from tag ") + _tag);
}

void sepRegFile::readComplexDoubleWindow(const std::vector<int> &nw,
                                   const std::vector<int> &fw,
                                   const std::vector<int> &jw,
                                   std::complex<double> *array) {
  std::shared_ptr<hypercube> hyper = getHyper();
  setDataType(DATA_COMPLEXDOUBLE);

  std::vector<int> ng = hyper->getNs();
  if (ng.size() > nw.size()) {
    for (int i = nw.size(); i < ng.size(); i++) {
      if (ng[i] > 1) error("number of dimension does not equal data size");
    }
  }
  if (nw.size() < ng.size() || fw.size() < ng.size() || jw.size() < jw.size()) {
    error("number of dimensions does not equal data size");
  }

  int ndim = ng.size();
  if (0 != sreed_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 16, array))
    error(std::string("trouble reading data from tag ") + _tag);
}


void sepRegFile::writeComplexStream(const std::complex<float> *array,
                                    const long long npts) {
  long long nptsT = npts * 8;
  long long ierr = srite_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble write from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}

void sepRegFile::writeComplexDoubleStream(const std::complex<double> *array,
                                    const long long npts) {
  set_format(_tag.c_str(), "native_double");
  long long nptsT = npts * 16;
  long long ierr = srite_big(_tag.c_str(), (void *)array, nptsT);
  if (ierr != nptsT)
    error(std::string("Trouble write from ") + _tag + std::string(" after ") +
          std::to_string(ierr) + std::string(" bytes"));
}


void sepRegFile::writeComplexWindow(const std::vector<int> &nw,
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
  int ndim = ng.size();
  if (0 != srite_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 8, array))
    error(std::string("trouble writing data to tag ") + _tag);
}
void sepRegFile::writeComplexDoubleWindow(const std::vector<int> &nw,
                                    const std::vector<int> &fw,
                                    const std::vector<int> &jw,
                                    const std::complex<double> *array) {
  setDataType(DATA_COMPLEXDOUBLE);

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
  int ndim = ng.size();
  if (0 != srite_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 16, array))
    error(std::string("trouble writing data to tag ") + _tag);
}
void sepRegFile::readByteWindow(const std::vector<int> &nw,
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
  int ndim = ng.size();
  if (0 != sreed_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 1, array))
    error(std::string("trouble reading data from tag ") + _tag);
}
void sepRegFile::writeFloatWindow(const std::vector<int> &nw,
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
  int ndim = ng.size();
  if (0 != srite_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 4, array))
    error(std::string("trouble writing data to tag ") + _tag);
}
void sepRegFile::writeByteWindow(const std::vector<int> &nw,
                                 const std::vector<int> &fw,
                                 const std::vector<int> &jw,
                                 const unsigned char *array) {
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
  int ndim = ng.size();
  if (0 != srite_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 1, array))
    error(std::string("trouble writing data to tag ") + _tag);
}

void sepRegFile::writeDoubleWindow(const std::vector<int> &nw,
                                   const std::vector<int> &fw,
                                   const std::vector<int> &jw,
                                   const double *array) {
  set_format(_tag.c_str(), "native_double");
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
  int ndim = ng.size();
  if (0 != srite_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 8, array))
    error(std::string("trouble writing data to tag ") + _tag);
}
void sepRegFile::writeIntWindow(const std::vector<int> &nw,
                                const std::vector<int> &fw,
                                const std::vector<int> &jw, const int *array) {
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
  int ndim = ng.size();
  if (0 != srite_window(_tag.c_str(), &ndim, ng.data(), nw.data(), fw.data(),
                        jw.data(), 4, array))
    error(std::string("trouble writing data to tag ") + _tag);
}

void sepRegFile::readDescription(const int ndimMax) {
  int ndim;
  int esize = getInt("esize", 4);
  if (esize == 1)
    setDataType(DATA_BYTE);
  else if (esize == 4) {
    std::string format =
        getString(std::string("data_format"), std::string("xdr_float"));
    if (format == std::string("xdr_float") ||
        format == std::string("native_float"))
      setDataType(DATA_FLOAT);
    else if (format == std::string("xdr_int") ||
             format == std::string("native_int"))
      setDataType(DATA_INT);
    else
      error(std::string("Unknown data type " + format));
  } else if (esize == 8) {
    std::string format =
        getString(std::string("data_format"), std::string("xdr_float"));
    if (format == std::string("xdr_float") || format == "native_float")
      setDataType(DATA_COMPLEX);
    else if (format == std::string("nativie_double"))
      setDataType(DATA_DOUBLE);
    else  // For now default to complex
      setDataType(DATA_COMPLEX);
} else if (esize == 16) {
    std::string format =std::string("nativie_double");
    setDataType(DATA_COMPLEXDOUBLE);
  } else
    error(std::string("Only know about esize=16,8, 4 or 1"));
  std::vector<axis> axes;
  sep_get_number_data_axes(_tag.c_str(), &ndim);
  if (ndimMax != -1 && ndimMax > ndim) ndim = ndimMax;
fprintf(stderr,"IN LOOP \n");
  putInt("esize",esize);
  for (int i = 1; i <= ndim; i++) {
    int n;
    float o, d;
    char label[1024];
    sep_get_data_axis_par(_tag.c_str(), &i, &n, &o, &d, label);
    axes.push_back(axis(n, o, d, std::string(label)));
  }
fprintf(stderr,"IN2LOOP \n");
  std::shared_ptr<hypercube> hyper(new hypercube(axes));
  setHyper(hyper);

}
void sepRegFile::writeDescription() {
  std::shared_ptr<hypercube> hyper = getHyper();
  std::vector<axis> axes = hyper->returnAxes(hyper->getNdim());

  for (int i = 1; i <= axes.size(); i++) {
    int n = axes[i - 1].n;
    float o = axes[i - 1].o;
    float d = axes[i - 1].d;
    char label[1024];
    std::copy(axes[i - 1].label.begin(), axes[i - 1].label.end(), label);
    label[axes[i - 1].label.length()] = '\0';
    sep_put_data_axis_par(_tag.c_str(), &i, &n, &o, &d, label);
  }
  for (int i = axes.size() + 1; i <= 8; i++) {
    int n = 1;
    float o = 0., d = 1.;
    sep_put_data_axis_par(_tag.c_str(), &i, &n, &o, &d, "none");
  }
  int esize = 4;

  switch (getDataType()) {
    case DATA_INT:
      set_format(_tag.c_str(), "xdr_int");
      break;
    case DATA_DOUBLE:
      set_format(_tag.c_str(), "native_float");
      esize = 8;
      break;
    case DATA_COMPLEX:
      set_format(_tag.c_str(), "xdr_int");
      esize = 8;
      break;
    case DATA_COMPLEXDOUBLE:
      set_format(_tag.c_str(), "native_double");
      esize = 16;
      break;
    case DATA_BYTE:
      set_format(_tag.c_str(), "xdr_byte");
      esize = 1;
      break;
    default:
      set_format(_tag.c_str(), "xdr_float");
      break;
  }

  auxputch("esize", "d", &esize, _tag.c_str());
}
void sepRegFile::message(const std::string &arg) const {
  sepwarn(0, arg.c_str());
}
void sepRegFile::seekTo(const long long iv, const int whence) {
  sseek(_tag.c_str(), iv, whence);
}
Json::Value sepRegFile::getDescription() {
  char *tmp_ch;
  int nsz = 100000;
  int nout = nsz * 2;
  tmp_ch = new char[1];
  grab_history(_tag.c_str(), tmp_ch, nsz, &nout);
  nsz = nout + 1;
  delete[] tmp_ch;
  tmp_ch = new char[nout];
  grab_history(_tag.c_str(), tmp_ch, nsz, &nout);
  std::string d = tmp_ch;
  delete[] tmp_ch;
  Json::Value j;
  j["SEPFILE"] = d;
  return j;
}

void sepRegFile::putDescription(const std::string &title,
                                const Json::Value &desc) {
  std::stringstream stream;
  stream << desc;
  std::string tmp = std::string("FROM ") + title;
  char delim = '\n';  // Ddefine the delimiter to split by
  auxputhead(_tag.c_str(), "%s\n", tmp.c_str());
  while (std::getline(stream, tmp, delim)) {
    auxputhead(_tag.c_str(), "%s\n", tmp.c_str());
  }
}
