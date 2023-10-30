#include <doubleHyper.h>
#include <hypercube.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tbb.h>
#include <iostream>
#include <random>
using namespace SEP;

void doubleHyper::add(const std::shared_ptr<doubleHyper> vec2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<doubleHyper> vec2H =
      std::dynamic_pointer_cast<doubleHyper>(vec2);
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] += vec2H->_vals[i];
                    });
  calcCheckSum();
}
void doubleHyper::mult(const std::shared_ptr<doubleHyper> vec2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<doubleHyper> vec2H =
      std::dynamic_pointer_cast<doubleHyper>(vec2);

  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] *= vec2H->_vals[i];
                    });
  calcCheckSum();
}
void doubleHyper::scaleAdd(std::shared_ptr<doubleHyper> vec2, const double sc1,
                           const double sc2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<doubleHyper> vec2H =
      std::dynamic_pointer_cast<doubleHyper>(vec2);
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] = _vals[i] * sc1 + sc2 * vec2H->_vals[i];
                    });
  calcCheckSum();
}
void doubleHyper::clipVector(const std::shared_ptr<doubleHyper> beg,
                             const std::shared_ptr<doubleHyper> end) {
  if (!checkSame(beg)) throw(std::string("Vectors not of the same space"));
  if (!checkSame(end)) throw(std::string("Vectors not of the same space"));

  std::shared_ptr<doubleHyper> begH =
      std::dynamic_pointer_cast<doubleHyper>(beg);
  std::shared_ptr<doubleHyper> endH =
      std::dynamic_pointer_cast<doubleHyper>(end);
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] = std::min(endH->_vals[i],
                                            std::max(begH->_vals[i], _vals[i]));
                    });
  calcCheckSum();
}
void doubleHyper::calcHisto(std::shared_ptr<int1DReg> &histo, float mn,
                            float mx) {
  long long nelem = histo->getHyper()->getN123();
  float delta = (mx - mn) / (nelem - 1);
  float idelta = 1. / delta;
  std::vector<int> h = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()),
      std::vector<int>(nelem, 0),
      [&](const tbb::blocked_range<size_t> &r, std::vector<int> tmp) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          int ielem = std::max(
              (long long)0,
              std::min(nelem - 1, (long long)((_vals[i] - mn) * idelta)));
          tmp[ielem] += 1;
        }
        return tmp;
      },
      [&](std::vector<int> tmp1, std::vector<int> tmp2) {
        std::vector<int> tmp = tmp1;
        for (int i = 0; i < nelem; i++) tmp[i] += tmp2[i];
        return tmp;
      });
  for (long long i = 0; i < nelem; i++) {
    histo->getVals()[i] = h[i];
  }
}
void doubleHyper::signum() {
  if (spaceOnly()) throw(std::string("Vectors not allocated"));
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i) {
                        if (_vals[i] > 1e-20)
                          _vals[i] = 1;
                        else if (_vals[i] < -1e-20)
                          _vals[i] = -1;
                        else
                          _vals[i] = 0;
                      }
                    });

  calcCheckSum();
}
double doubleHyper::cent(const long long iv, const int js) const {
  long long n = getHyper()->getN123() / js;
  double *x = new double[n];
  const double *in = getCVals();
  if (js < 1)
    throw SEPException(std::string("j must be positive ") + std::to_string(js));
  for (auto i = 0; i < n; i++) {
    x[i] = in[i * js];
  }
  double w = (double)(iv) / (double)n;
  long long q =
      std::max((long long)0, std::min((long long)(n - 1), (long long)(w * n)));
  double *i, *j, ak;
  double *low, *hi, buf, *k;
  for (low = x, hi = x + n - 1, k = x + q; low < hi;) {
    ak = *k;
    i = low;
    j = hi;
    do {
      while (*i < ak) i++;
      while (*j > ak) j--;
      if (i <= j) {
        buf = *i;
        *i++ = *j;
        *j-- = buf;
      }
    } while (i <= j);
    if (j < k) low = i;
    if (k < i) hi = j;
  }
  double vv = *k;
  delete[] x;
  return vv;
}
void doubleHyper::clip(const double bclip, const double eclip, bool outer) {
  if (spaceOnly()) throw(std::string("Vectors not allocated"));
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      if (outer) {
                        for (long long i = r.begin(); i != r.end(); ++i)
                          _vals[i] = std::min(eclip, std::max(bclip, _vals[i]));
                      } else {
                        for (long long i = r.begin(); i != r.end(); ++i) {
                          if (_vals[i] < eclip && _vals[i] > bclip)
                            _vals[i] = 0;
                        }
                      }
                    });
}
void doubleHyper::scale(double sc) {
  if (spaceOnly()) throw(std::string("Vectors not allocated"));
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] = _vals[i] * sc;
                    });
  calcCheckSum();
}
void doubleHyper::random() {
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] = ((double)rand() / (RAND_MAX)) - .5;
                    });
  calcCheckSum();
}

double doubleHyper::norm(const int n) const {
  double dt = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), 0.,
      [&](const tbb::blocked_range<size_t> &r, double v) {
        if (n == 1) {
          for (size_t i = r.begin(); i != r.end(); ++i) {
            v += (double)fabs(_vals[i]);
          }
        } else {
          for (size_t i = r.begin(); i != r.end(); ++i) {
            v += (double)_vals[i] * (double)_vals[i];
          }
        }
        return v;
      },
      [](double a, double b) { return a + b; });
  if (n == 2) return sqrtf(dt);
  return dt;
}
void doubleHyper::set(const double val) {
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] = val;
                    });
  calcCheckSum();
}
double doubleHyper::dot(const std::shared_ptr<doubleHyper> vec2) const {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<doubleHyper> vec2H =
      std::dynamic_pointer_cast<doubleHyper>(vec2);

  double dot = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), 0.,
      [&](const tbb::blocked_range<size_t> &r, double v) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          v += (double)_vals[i] * (double)vec2H->_vals[i];
        }
        return v;
      },
      [](double a, double b) { return a + b; });

  return dot;
}
void doubleHyper::createMask(const float zero, const float err) {
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i) {
                        if (fabs(_vals[i] - zero) > err)

                          _vals[i] = 0.;
                        else
                          _vals[i] = 1;
                      }
                    });
  calcCheckSum();
}

void doubleHyper::infoStream(const int lev, std::stringstream &x) {
  getHyper()->infoStream(x);
  if (spaceOnly())
    x << "Only space\n";
  else {
    x << "Allocated\n";
    long long npts = std::min((const long long)lev, getHyper()->getN123());
    for (long long i = 0; i < npts; i++)
      x << std::to_string(i) << std::string(" ") << std::to_string(_vals[i])
        << std::endl;
  }
}
void doubleHyper::softClip(const float scale) {
  float sc2 = scale * scale;
  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] = scale * _vals[i] /
                                   sqrtf(1. + sc2 * _vals[i] * _vals[i]);
                    });
  calcCheckSum();
}

double doubleHyper::absMax() const {
  double val = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), -1.e31,
      [&](const tbb::blocked_range<size_t> &r, double v) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          v = std::max(v, fabs(_vals[i]));
        }
        return v;
      },
      [](double a, double b) {
        if (a > b) return a;
        return b;
      });
  return (double)val;
}
double doubleHyper::max() const {
  double val = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), -1.e31,
      [&](const tbb::blocked_range<size_t> &r, double v) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          v = std::max(v, _vals[i]);
        }
        return v;
      },
      [](double a, double b) {
        if (a > b) return a;
        return b;
      });
  return (double)val;
}
double doubleHyper::min() const {
  double val = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), 1.e31,
      [&](const tbb::blocked_range<size_t> &r, double v) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          v = std::min(v, _vals[i]);
        }
        return v;
      },
      [](double a, double b) {
        if (a > b) return b;
        return a;
      });
  return (double)val;
}
void doubleHyper::calcCheckSum() {
  uint32_t sum1 = 0, sum2 = 0;
  uint32_t *data = (uint32_t *)_vals;
  uint32_t mx = 4294967295;
  for (long long i = 0; i < getHyper()->getN123(); i++) {
    sum1 = (sum1 + data[i]) % mx;
    sum2 = (sum2 + sum1) % mx;
  }
  setCheckSum(sum2 * 2 ^ 32 + sum1);
}

bool doubleHyper::checkSame(const std::shared_ptr<doubleHyper> vec2) const {
  if (!vec2) {
    throw SEPException("Vec2 is not allocated");
    return false;
  }
  bool b;
  try {
    b = getHyper()->checkSame(vec2->getHyper());
  } catch (SEPException &e) {
    throw SEPException(e.getMessage());
  }

  // return b;
  return true;
}
