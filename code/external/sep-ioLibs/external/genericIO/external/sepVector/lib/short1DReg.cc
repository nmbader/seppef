#include <short1DReg.h>
#include "SEPException.h"
using namespace SEP;
std::shared_ptr<short1DReg> short1DReg::clone() const {
  if (getSpaceOnly()) {
    std::shared_ptr<short1DReg> x(new short1DReg(getHyper()));

    return x;
  } else {
    std::shared_ptr<short1DReg> x(new short1DReg(getHyper(), *_mat));
    return x;
  }
}
std::shared_ptr<short1DReg> short1DReg::cloneSpace() const {
  std::shared_ptr<short1DReg> x(new short1DReg(getHyper()));
  x->_mat = 0;
  x->setSpace();

  return x;
}
void short1DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  setHyper(hyp);

  _vecType = "vec 1d int";
  const std::vector<SEP::axis> axes = hyp->getAxes();
  if (axes.size() != 1) throw(SEPException("Hypercube must by 1-D"));
  _mat.reset(new short1D(boost::extents[axes[0].n]));
  setData(_mat->data());
}
void short1DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                          const short1D &vals) {
  setHyper(hyp);

  const std::vector<SEP::axis> axes = hyp->getAxes();
  if (axes.size() != 1) throw(SEPException("Hypercube must by 1-D"));
  if (axes[0].n != vals.shape()[0])
    throw(SEPException(std::string("SHort array=") +
                       std::to_string(vals.shape()[0]) +
                       std::string(" must be the same size as data=") +
                       std::to_string(axes[0].n)));

  _mat.reset(new short1D(boost::extents[axes[0].n]));
  setData(_mat->data());
  for (long long i = 0; i < axes[0].n; i++) (*_mat)[i] = vals[i];
}
