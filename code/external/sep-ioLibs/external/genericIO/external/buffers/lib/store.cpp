#include <assert.h>
#include <store.h>
#include <iostream>
#include "SEPException.h"
using namespace SEP::IO;

storeInt::storeInt(const size_t n, void *buf) {
  _buf = new int[n];
  _n = n;
  memcpy(_buf, buf, n * sizeof(int));
}
void storeInt::getData(std::shared_ptr<storeBase> buf) const {
  std::shared_ptr<storeInt> b = std::dynamic_pointer_cast<storeInt>(buf);
  if (!b) throwError("storeInt", returnStorageType(buf));
  assert(b->getSize() <= getSize());
  memcpy(b->getPtr(), _buf, getSize() * sizeof(int));
}

void storeInt::putData(const std::shared_ptr<storeBase> buf) {
  const std::shared_ptr<storeInt> b = std::dynamic_pointer_cast<storeInt>(buf);
  if (!b) throwError("storeInt", returnStorageType(buf));
  assert(b->getSize() <= getSize());
  memcpy(_buf, b->getPtr(), getSize() * sizeof(int));
}
std::shared_ptr<storeBase> storeInt::clone() const {
  std::shared_ptr<storeInt> m(new storeInt((int)getSize()));
  memcpy(m->_buf, _buf, getSize() * sizeof(int));

  return m;
}
void storeInt::getWindow(const std::vector<int> &nwL,
                         const std::vector<int> &fwL,
                         const std::vector<int> &jwL,
                         const std::vector<int> &nbL,
                         const std::vector<int> &fwG,
                         const std::vector<int> &nbG, void *bufIn) {
  int *buf = (int *)bufIn;

  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                buf[i0L + f1G] = _buf[f1L + i0L * jwL[0]];
              }
            }
          }
        }
      }
    }
  }
}
void storeInt::putWindow(const std::vector<int> &nwL,
                         const std::vector<int> &fwL,
                         const std::vector<int> &jwL,
                         const std::vector<int> &nbL,
                         const std::vector<int> &fwG,
                         const std::vector<int> &nbG, const void *bufIn) {
  const int *buf = (int *)bufIn;
  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                _buf[f1L + i0L * jwL[0]] = buf[i0L + f1G];
              }
            }
          }
        }
      }
    }
  }
}

storeByte::storeByte(const size_t n, void *buf) {
  _buf = new unsigned char[n];
  _n = n;

  memcpy(_buf, buf, n * sizeof(unsigned char));
}
void storeByte::getData(std::shared_ptr<storeBase> buf) const {
  std::shared_ptr<storeByte> b = std::dynamic_pointer_cast<storeByte>(buf);
  if (!b) throwError("storeByte", returnStorageType(buf));
  assert(b->getSize() <= getSize());
  memcpy(b->_buf, _buf, getSize() * sizeof(unsigned char));
}

void storeByte::putData(const std::shared_ptr<storeBase> buf) {
  const std::shared_ptr<storeByte> b =
      std::dynamic_pointer_cast<storeByte>(buf);
  if (!b) throwError("storeByte", returnStorageType(buf));
  assert(b->getSize() == getSize());
  memcpy(_buf, b->_buf, getSize() * sizeof(unsigned char));
}
std::shared_ptr<storeBase> storeByte::clone() const {
  std::shared_ptr<storeByte> m(new storeByte((int)getSize()));
  memcpy(m->_buf, _buf, getSize() * sizeof(unsigned char));

  return m;
}
void storeByte::getWindow(const std::vector<int> &nwL,
                          const std::vector<int> &fwL,
                          const std::vector<int> &jwL,
                          const std::vector<int> &nbL,
                          const std::vector<int> &fwG,
                          const std::vector<int> &nbG, void *bufIn) {
  unsigned char *buf = (unsigned char *)bufIn;
  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                buf[i0L + f1G] = _buf[f1L + i0L * jwL[0]];
              }
            }
          }
        }
      }
    }
  }
}
void storeByte::putWindow(const std::vector<int> &nwL,
                          const std::vector<int> &fwL,
                          const std::vector<int> &jwL,
                          const std::vector<int> &nbL,
                          const std::vector<int> &fwG,
                          const std::vector<int> &nbG, const void *bufIn) {
  const unsigned char *buf = (unsigned char *)bufIn;
  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                _buf[f1L + i0L * jwL[0]] = buf[i0L + f1G];
              }
            }
          }
        }
      }
    }
  }
}

storeFloat::storeFloat(const size_t n, void *buf) {
  _buf = new float[n];
  _n = n;

  memcpy(_buf, buf, n * sizeof(float));
}
void storeFloat::getData(std::shared_ptr<storeBase> buf) const {
  std::shared_ptr<storeFloat> b = std::dynamic_pointer_cast<storeFloat>(buf);
  if (!b) throwError("storeFloat", returnStorageType(buf));
  assert(b->getSize() <= getSize());
  memcpy(b->_buf, _buf, getSize() * sizeof(float));
}

void storeFloat::putData(const std::shared_ptr<storeBase> buf) {
  const std::shared_ptr<storeFloat> b =
      std::dynamic_pointer_cast<storeFloat>(buf);
  if (!b) throwError("storeFloat", returnStorageType(buf));
  assert(b->getSize() == getSize());
  memcpy(_buf, b->_buf, getSize() * sizeof(float));
}
std::shared_ptr<storeBase> storeFloat::clone() const {
  std::shared_ptr<storeFloat> m(new storeFloat((int)getSize()));
  memcpy(m->_buf, _buf, getSize() * sizeof(float));

  return m;
}
void storeFloat::getWindow(const std::vector<int> &nwL,
                           const std::vector<int> &fwL,
                           const std::vector<int> &jwL,
                           const std::vector<int> &nbL,
                           const std::vector<int> &fwG,

                           const std::vector<int> &nbG, void *bufIn) {
  float *buf = (float *)bufIn;

  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
	      //if(jwL[0]==1){
	//	       memcpy(&buf[f1G],&buf[f1L],nwL[0]*sizeof(float));
	 //     }
	  //    else{
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                buf[i0L + f1G] = _buf[f1L + i0L * jwL[0]];
              }
	   //   }
            }
          }
        }
      }
    }
  }
}
void storeFloat::info(const std::string &v) const {
  double mn = 1e31, mx = -1e31, sm = 0;
  for (auto i = 0; i < getSize(); i++) {
    mn = std::min((float)mn, _buf[i]);
    mx = std::max((float)mx, _buf[i]);
    sm += fabs(_buf[i]);
  }
  std::cerr << v << " size=" << getSize() << " min=" << mn << " max=" << mx
            << " sm=" << sm << std::endl;
}
void storeFloat::putWindow(const std::vector<int> &nwL,
                           const std::vector<int> &fwL,
                           const std::vector<int> &jwL,
                           const std::vector<int> &nbL,
                           const std::vector<int> &fwG,
                           const std::vector<int> &nbG, const void *bufIn) {

  const float *buf = (float *)bufIn;
  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
	      //if(jwL[0]==1){
	//	       memcpy(&_buf[f1L],&buf[f1G],nwL[0]*sizeof(float));
	 //     }
	  //    else{
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                _buf[f1L + i0L * jwL[0]] = buf[i0L + f1G];
           //   }
	      }
            }
          }
        }
      }
    }
  }
}

storeDouble::storeDouble(const size_t n, void *buf) {
  _buf = new double[n];
  _n = n;

  memcpy(_buf, buf, n * sizeof(float));
}
void storeDouble::getData(std::shared_ptr<storeBase> buf) const {
  std::shared_ptr<storeDouble> b = std::dynamic_pointer_cast<storeDouble>(buf);
  if (!b) throwError("storeDouble", returnStorageType(buf));
  assert(b->getSize() <= getSize());
  memcpy(b->_buf, _buf, getSize() * sizeof(double));
}

void storeDouble::putData(const std::shared_ptr<storeBase> buf) {
  const std::shared_ptr<storeDouble> b =
      std::dynamic_pointer_cast<storeDouble>(buf);
  if (!b) throwError("storeDouble", returnStorageType(buf));
  assert(b->getSize() == getSize());
  memcpy(_buf, b->_buf, getSize() * sizeof(double));
}
std::shared_ptr<storeBase> storeDouble::clone() const {
  std::shared_ptr<storeDouble> m(new storeDouble((int)getSize()));
  memcpy(m->_buf, _buf, getSize() * sizeof(double));

  return m;
}
void storeDouble::getWindow(const std::vector<int> &nwL,
                            const std::vector<int> &fwL,
                            const std::vector<int> &jwL,
                            const std::vector<int> &nbL,
                            const std::vector<int> &fwG,
                            const std::vector<int> &nbG, void *bufIn) {
  double *buf = (double *)bufIn;

  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                buf[i0L + f1G] = _buf[f1L + i0L * jwL[0]];
              }
            }
          }
        }
      }
    }
  }
}
void storeDouble::putWindow(const std::vector<int> &nwL,
                            const std::vector<int> &fwL,
                            const std::vector<int> &jwL,
                            const std::vector<int> &nbL,
                            const std::vector<int> &fwG,
                            const std::vector<int> &nbG, const void *bufIn) {
  const double *buf = (double *)bufIn;
  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                _buf[f1L + i0L * jwL[0]] = buf[i0L + f1G];
              }
            }
          }
        }
      }
    }
  }
}

storeComplex::storeComplex(const size_t n, void *buf) {
  _buf = new std::complex<float>[n];
  _n = n;

  memcpy(_buf, buf, n * 2 * sizeof(float));
}
void storeComplex::getData(std::shared_ptr<storeBase> buf) const {
  std::shared_ptr<storeComplex> b =
      std::dynamic_pointer_cast<storeComplex>(buf);
  if (!b) throwError("storeComplex", returnStorageType(buf));

  assert(b);
  assert(b->getSize() <= getSize());
  memcpy(b->_buf, _buf, getSize() * sizeof(std::complex<float>));
}

void storeComplex::putData(const std::shared_ptr<storeBase> buf) {
  const std::shared_ptr<storeComplex> b =
      std::dynamic_pointer_cast<storeComplex>(buf);
  if (!b) throwError("storeComplex", returnStorageType(buf));
  assert(b->getSize() == getSize());
  memcpy(_buf, b->_buf, getSize() * sizeof(std::complex<float>));
}
std::shared_ptr<storeBase> storeComplex::clone() const {
  std::shared_ptr<storeComplex> m(new storeComplex((int)getSize()));
  memcpy(m->_buf, _buf, getSize() * sizeof(std::complex<float>));

  return m;
}
void storeComplex::getWindow(const std::vector<int> &nwL,
                             const std::vector<int> &fwL,
                             const std::vector<int> &jwL,
                             const std::vector<int> &nbL,
                             const std::vector<int> &fwG,
                             const std::vector<int> &nbG, void *bufIn) {
  std::complex<float> *buf = (std::complex<float> *)bufIn;

  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                buf[i0L + f1G] = _buf[f1L + i0L * jwL[0]];
              }
            }
          }
        }
      }
    }
  }
}
void storeComplex::putWindow(const std::vector<int> &nwL,
                             const std::vector<int> &fwL,
                             const std::vector<int> &jwL,
                             const std::vector<int> &nbL,
                             const std::vector<int> &fwG,
                             const std::vector<int> &nbG, const void *bufIn) {
  const std::complex<float> *buf = (std::complex<float> *)bufIn;
  for (int i6L = 0; i6L < nwL[6]; i6L++) {
    size_t f6L = nbL[6] * (fwL[6] + i6L * jwL[6]);
    size_t f6G = nbG[6] * (fwG[6] + i6L);
    for (int i5L = 0; i5L < nwL[5]; i5L++) {
      size_t f5L = f6L + nbL[5] * (fwL[5] + i5L * jwL[5]);
      size_t f5G = nbG[5] * (fwG[5] + i5L) + f6G;
      for (int i4L = 0; i4L < nwL[4]; i4L++) {
        size_t f4L = f5L + nbL[4] * (fwL[4] + jwL[4] * i4L);
        size_t f4G = nbG[4] * (fwG[4] + i4L) + f5G;
        for (int i3L = 0; i3L < nwL[3]; i3L++) {
          size_t f3L = f4L + nbL[3] * (fwL[3] + jwL[3] * i3L);
          size_t f3G = nbG[3] * (fwG[3] + i3L) + f4G;
          for (int i2L = 0; i2L < nwL[2]; i2L++) {
            size_t f2L = f3L + nbL[2] * (fwL[2] + jwL[2] * i2L);
            size_t f2G = nbG[2] * (fwG[2] + i2L) + f3G;
            for (int i1L = 0; i1L < nwL[1]; i1L++) {
              size_t f1L = f2L + nbL[1] * (fwL[1] + jwL[1] * i1L) + fwL[0];
              size_t f1G = nbG[1] * (fwG[1] + i1L) + f2G + fwG[0];
              for (size_t i0L = 0; i0L < nwL[0]; i0L++) {
                _buf[f1L + i0L * jwL[0]] = buf[i0L + f1G];
              }
            }
          }
        }
      }
    }
  }
}

std::shared_ptr<storeBase> SEP::IO::returnStorage(const dataType state,
                                                  const size_t n) {
  switch (state) {
    case DATA_INT: {
      std::shared_ptr<storeInt> x(new storeInt(n));
      return x;
    } break;
    case DATA_FLOAT: {
      std::shared_ptr<storeFloat> y(new storeFloat(n));
      return y;
    } break;
    case DATA_COMPLEX: {
      std::shared_ptr<storeComplex> z(new storeComplex(n));
      return z;
    } break;
    case DATA_DOUBLE: {
      std::shared_ptr<storeDouble> a(new storeDouble(n));
      return a;
    } break;
    case DATA_BYTE: {
      std::shared_ptr<storeByte> b(new storeByte(n));
      return b;
    } break;
    default:
      throw SEPException("Unknown data type");
  }
}

std::string SEP::IO::returnStorageType(const std::shared_ptr<storeBase> bufIn) {
  {
    const std::shared_ptr<storeByte> by =
        std::dynamic_pointer_cast<storeByte>(bufIn);
    if (by) return "storeByte";
  }
  {
    const std::shared_ptr<storeInt> by =
        std::dynamic_pointer_cast<storeInt>(bufIn);
    if (by) return "storeInt";
  }
  {
    const std::shared_ptr<storeDouble> by =
        std::dynamic_pointer_cast<storeDouble>(bufIn);

    if (by) return "storeDouble";
  }
  {
    const std::shared_ptr<storeFloat> by =
        std::dynamic_pointer_cast<storeFloat>(bufIn);

    if (by) return "storeFloat";
  }
  {
    const std::shared_ptr<storeComplex> by =
        std::dynamic_pointer_cast<storeComplex>(bufIn);
    if (by) return "storeComplex";
  }
  return "storeUnknown";
}

void SEP::IO::throwError(const std::string myType, const std::string typeSent) {
  throw SEPException(std::string("Wrong type passed Expected:") + myType +
                     std::string(" got ") + std::string(typeSent));
}
