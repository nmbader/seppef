#include <complexDoubleHyper.h>
#include <hypercube.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tbb.h>
#include <iostream>
#include <random>

using namespace SEP;


void complexDoubleHyper::add(const std::shared_ptr<complexDoubleHyper> vec2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<complexDoubleHyper> vec2H =
      std::dynamic_pointer_cast<complexDoubleHyper>(vec2);

  for (long long i = 0; i < getHyper()->getN123(); i++)

    tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                      [&](const tbb::blocked_range<long long> &r) {
                        for (long long i = r.begin(); i != r.end(); ++i)
                          _vals[i] += vec2H->_vals[i];
                      });
  calcCheckSum();
}

void complexDoubleHyper::scaleAdd(std::shared_ptr<complexDoubleHyper> vec2,
                            const double sc1, const double sc2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));

  std::shared_ptr<complexDoubleHyper> vec2H =
      std::dynamic_pointer_cast<complexDoubleHyper>(vec2);

  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i) {
                        std::complex<double> f1 = {_vals[i].real() * (double)sc1,
                                                  _vals[i].imag() * (double)sc1};
                        std::complex<double> f2 = {
                            (double)sc2 * vec2H->_vals[i].real(),
                            (double)sc2 * vec2->_vals[i].imag()};
                        _vals[i] = f1 + f2;
                      };
                    });
  calcCheckSum();
}
void complexDoubleHyper::mult(const std::shared_ptr<complexDoubleHyper> vec2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<complexDoubleHyper> vec2H =
      std::dynamic_pointer_cast<complexDoubleHyper>(vec2);

  tbb::parallel_for(tbb::blocked_range<long long>(0, getHyper()->getN123()),
                    [&](const tbb::blocked_range<long long> &r) {
                      for (long long i = r.begin(); i != r.end(); ++i)
                        _vals[i] *= vec2H->_vals[i];
                    });
  calcCheckSum();
}
std::complex<double> complexDoubleHyper::dot(const std::shared_ptr<complexDoubleHyper> vec2) const {
  if (!checkSame(vec2)) throw(std::string("Vectors not of the same space"));
  std::shared_ptr<complexDoubleHyper> vec2H =
      std::dynamic_pointer_cast<complexDoubleHyper>(vec2);

  double dot = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, getHyper()->getN123()), 0.,
      [&](const tbb::blocked_range<size_t> &r, double v) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          v += (double)_vals[i].real() * (double)vec2H->_vals[i].real() +
               (double)_vals[i].imag() * (double)vec2H->_vals[i].imag();
        }
        return v;
      },
      [](double a, double b) { return a + b; });

  return dot;
}



void complexDoubleHyper::clipVector(const std::shared_ptr<doubleHyper> begV,
                              const std::shared_ptr<doubleHyper> endV) {
  if (!getHyper()->checkSame(begV->getHyper()))
    throw(std::string("Vectors not of the same space"));
  if (!getHyper()->checkSame(endV->getHyper()))
    throw(std::string("Vectors not of the same space"));

  std::shared_ptr<doubleHyper> begH =
      std::dynamic_pointer_cast<doubleHyper>(begV);
  std::shared_ptr<doubleHyper> endH =
      std::dynamic_pointer_cast<doubleHyper>(endV);

  double *beg = begH->getVals(), *end = endH->getVals();
  tbb::parallel_for(
      tbb::blocked_range<long long>(0, getHyper()->getN123()),
      [&](const tbb::blocked_range<long long> &r) {
        for (long long i = r.begin(); i != r.end(); ++i)
          _vals[i] =
              std::min(end[i], std::max(beg[i], (double)sqrtf(std::norm(_vals[i]))));
      });
  calcCheckSum();
}

void complexDoubleHyper::clipVector(const std::shared_ptr<complexDoubleHyper> beg,
                              const std::shared_ptr<complexDoubleHyper> end) {
  if (!checkSame(beg)) throw(std::string("Vectors not of the same space"));
  if (!checkSame(end)) throw(std::string("Vectors not of the same space"));

  std::shared_ptr<complexDoubleHyper> begH =
      std::dynamic_pointer_cast<complexDoubleHyper>(beg);
  std::shared_ptr<complexDoubleHyper> endH =
      std::dynamic_pointer_cast<complexDoubleHyper>(end);
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

void complexDoubleHyper::scale(const double sc) {
  tbb::parallel_for(
      tbb::blocked_range<long long>(0, getHyper()->getN123()),
      [&](const tbb::blocked_range<long long> &r) {
        for (long long i = r.begin(); i != r.end(); ++i)
          _vals[i] = {_vals[i].real() * (double)sc, _vals[i].imag() * (double)sc};
      });
  calcCheckSum();
}

void complexDoubleHyper::random() {
  if (getSpaceOnly()) throw(std::string("Vectors not allocated"));
  for (long long ii = 0; ii < getHyper()->getN123(); ii++)
    _vals[ii] = {(double)((double)rand() / (RAND_MAX)-.5),
                 (double)((double)rand() / (RAND_MAX)-.5)};
  calcCheckSum();
}

double complexDoubleHyper::norm(const int n) const {
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

void complexDoubleHyper::set(const std::complex<double> val) {
  for (long long i = 0; i < getHyper()->getN123(); i++) _vals[i] = val;
  calcCheckSum();
}

void complexDoubleHyper::infoStream(const int lev, std::stringstream &x) {
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

void complexDoubleHyper::calcCheckSum() {
  uint32_t sum1 = 0, sum2 = 0;
  uint32_t *data = (uint32_t *)_vals;
  uint32_t mx = 4294967295;
  for (long long i = 0; i < getHyper()->getN123(); i++) {
    sum1 = (sum1 + data[i]) % mx;
    sum2 = (sum2 + sum1) % mx;
  }
  setCheckSum(sum2 * 2 ^ 32 + sum1);
}

bool complexDoubleHyper::checkSame(const std::shared_ptr<complexDoubleHyper> vec2) const {
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
