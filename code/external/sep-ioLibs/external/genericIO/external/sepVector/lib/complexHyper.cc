#include <complexHyper.h>
#include <hypercube.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tbb.h>
#include <iostream>
#include <random>
using namespace SEP;


void complexHyper::add(const std::shared_ptr<complexHyper> vec2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<complexHyper> vec2H =
      std::dynamic_pointer_cast<complexHyper>(vec2);

  for (long long i = 0; i < getHyper()->getN123(); i++)

    tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                      [&](const tbb::blocked_range<long long> &r) {
                        for (long long i = r.begin(); i != r.end(); ++i)
                          _vals[i] += vec2H->_vals[i];
                      });
  calcCheckSum();
}

void complexHyper::scaleAdd(std::shared_ptr<complexHyper> vec2,
                            const double sc1, const double sc2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));

  std::shared_ptr<complexHyper> vec2H =
      std::dynamic_pointer_cast<complexHyper>(vec2);

  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i) {
                        std::complex<float> f1 = {_vals[i].real() * (float)sc1,
                                                  _vals[i].imag() * (float)sc1};
                        std::complex<float> f2 = {
                            (float)sc2 * vec2H->_vals[i].real(),
                            (float)sc2 * vec2->_vals[i].imag()};
                        _vals[i] = f1 + f2;
                      };
                    });
  calcCheckSum();
}
void complexHyper::mult(const std::shared_ptr<complexHyper> vec2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<complexHyper> vec2H =
      std::dynamic_pointer_cast<complexHyper>(vec2);

  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] *= vec2H->_vals[i];
                    });
  calcCheckSum();
}
std::complex<double> complexHyper::dot(const std::shared_ptr<complexHyper> vec2) const {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<complexHyper> vec2H =
      std::dynamic_pointer_cast<complexHyper>(vec2);

  std::complex<double> dt = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), std::complex<double>(0.,0.),
      [&](const tbb::blocked_range<size_t> &r, std::complex<double> v) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
           v+=std::complex<double>( (double)_vals[i].real() * (double)vec2H->_vals[i].real() +
               (double)_vals[i].imag() * (double)vec2H->_vals[i].imag(),
               (double)_vals[i].real() * (double)vec2H->_vals[i].imag() -
               (double)_vals[i].imag() * (double)vec2H->_vals[i].real()
               );
        }
        return v;
      },
      [](std::complex<double> a, std::complex<double> b) { return a + b; });

  return dt;
}
  


void complexHyper::clipVector(const std::shared_ptr<floatHyper> begV,
                              const std::shared_ptr<floatHyper> endV) {
  if (!getHyper()->checkSame(begV->getHyper()))
    throw(std::string("Vectors not of the same space"));
  if (!getHyper()->checkSame(endV->getHyper()))
    throw(std::string("Vectors not of the same space"));

  std::shared_ptr<floatHyper> begH =
      std::dynamic_pointer_cast<floatHyper>(begV);
  std::shared_ptr<floatHyper> endH =
      std::dynamic_pointer_cast<floatHyper>(endV);

  float *beg = begH->getVals(), *end = endH->getVals();
  tbb::parallel_for(
      tbb::blocked_range<long long>(0, getHyper()->getN123()),
      [&](const tbb::blocked_range<long long> &r) {
        for (long long i = r.begin(); i != r.end(); ++i)
          _vals[i] =
              std::min(end[i], std::max(beg[i], sqrtf(std::norm(_vals[i]))));
      });
  calcCheckSum();
}

void complexHyper::clipVector(const std::shared_ptr<complexHyper> beg,
                              const std::shared_ptr<complexHyper> end) {
  if (!checkSame(beg)) throw(std::string("Vectors not of the same space"));
  if (!checkSame(end)) throw(std::string("Vectors not of the same space"));

  std::shared_ptr<complexHyper> begH =
      std::dynamic_pointer_cast<complexHyper>(beg);
  std::shared_ptr<complexHyper> endH =
      std::dynamic_pointer_cast<complexHyper>(end);
  tbb::parallel_for(
      tbb::blocked_range<long long>(0, getHyper()->getN123()),
      [&](const tbb::blocked_range<long long> &r) {
        for (long long i = r.begin(); i != r.end(); ++i)
          _vals[i] = {
              std::min(endH->_vals[i].real(),
                       std::max(begH->_vals[i].real(), _vals[i].real())),
              std::min(endH->_vals[i].imag(),
                       std::max(begH->_vals[i].imag(), _vals[i].imag()))

          };
      });
  calcCheckSum();
}

void complexHyper::scale(const double sc) {
  tbb::parallel_for(
      tbb::blocked_range<long long>(0, getHyper()->getN123()),
      [&](const tbb::blocked_range<long long> &r) {
        for (long long i = r.begin(); i != r.end(); ++i)
          _vals[i] = {_vals[i].real() * (float)sc, _vals[i].imag() * (float)sc};
      });
  calcCheckSum();
}

void complexHyper::random() {
  if (getSpaceOnly()) throw(std::string("Vectors not allocated"));
  for (long long ii = 0; ii < getHyper()->getN123(); ii++)
    _vals[ii] = {(float)((double)rand() / (RAND_MAX)-.5),
                 (float)((double)rand() / (RAND_MAX)-.5)};
  calcCheckSum();
}

double complexHyper::norm(const int n) const {
  double dt = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), 0.,
      [&](const tbb::blocked_range<size_t> &r, double v) {
        if (n == 1) {
          for (size_t i = r.begin(); i != r.end(); ++i) {
            v += (double)sqrt(std::norm(_vals[i]));
          }
        } else {
          for (size_t i = r.begin(); i != r.end(); ++i) {
            v += (double)std::norm(_vals[i]);
          }
        }
        return v;
      },
      [](double a, double b) { return a + b; });
  if (n == 2) return sqrtf(dt);

  return dt;
}

void complexHyper::set(const std::complex<float> val) {
  for (long long i = 0; i < getHyper()->getN123(); i++) _vals[i] = val;
  calcCheckSum();
}

void complexHyper::infoStream(const int lev, std::stringstream &x) {
  getHyper()->infoStream(x);
  if (getSpaceOnly())
    x << "Only space\n";
  else {
    x << "Allocated\n";
    long long npts = std::min((const long long)lev, getHyper()->getN123());
    //  for (long long i = 0; i < npts; i++)
    //  x << std::to_string(i) << std::string(" ") <<
    //  std::to_string(_vals[i])
    //  << std::endl;
  }
}

void complexHyper::calcCheckSum() {
  uint32_t sum1 = 0, sum2 = 0;
  uint32_t *data = (uint32_t *)_vals;
  uint32_t mx = 4294967295;
  for (long long i = 0; i < getHyper()->getN123(); i++) {
    sum1 = (sum1 + data[i]) % mx;
    sum2 = (sum2 + sum1) % mx;
  }
  setCheckSum(sum2 * 2 ^ 32 + sum1);
}

bool complexHyper::checkSame(const std::shared_ptr<complexHyper> vec2) const {
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

  return b;
  return true;
}
