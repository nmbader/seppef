#include <complex3DReg.h>
using namespace SEP;

std::shared_ptr<complex3DReg> complex3DReg::clone() const {
  if (getSpaceOnly()) {
    std::shared_ptr<complex3DReg> x(new complex3DReg(getHyper()));
    return x;
  } else {
    std::shared_ptr<complex3DReg> x(new complex3DReg(getHyper(), *_mat));

    return x;
  }
}
std::shared_ptr<complex3DReg> complex3DReg::cloneSpace() const {
  std::shared_ptr<complex3DReg> x(new complex3DReg(getHyper()));
  x->_mat = 0;
  x->setSpace();

  return x;
}

void complex3DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  if (axes.size() != 3)
    throw(SEPException(std::string("Axes size must be 3 is ") +
                       std::to_string(axes.size())));

  _mat.reset(new complex3D(boost::extents[axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
}
void complex3DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                            const complex3D &vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);
  if (axes.size() != 3)
    throw(SEPException(std::string("Axes size must be 3 is ") +
                       std::to_string(axes.size())));

  if (axes[0].n != vals.shape()[2])
    throw(SEPException(std::string("Axis 1 not the same (") +
                       std::to_string(axes[0].n) + std::string(",") +
                       std::to_string(vals.shape()[2]) + std::string(")")));
  if (axes[1].n != vals.shape()[1])
    throw(SEPException(std::string("Axis 2 not the same (") +
                       std::to_string(axes[1].n) + std::string(",") +
                       std::to_string(vals.shape()[1]) + std::string(")")));
  if (axes[2].n != vals.shape()[0])
    throw(SEPException(std::string("Axis 3 not the same (") +
                       std::to_string(axes[2].n) + std::string(",") +
                       std::to_string(vals.shape()[0]) + std::string(")")));
  _mat.reset(new complex3D(boost::extents[axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (long long k = 0; k < axes[2].n; k++) {
    for (long long j = 0; j < axes[1].n; j++) {
      for (long long i = 0; i < axes[0].n; i++) {
        (*_mat)[k][j][i] = vals[k][j][i];
      }
    }
  }
}
void complex3DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                            std::shared_ptr<complex3D> vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  if (axes.size() != 3)
    throw(SEPException(std::string("Axes size must be 3 is ") +
                       std::to_string(axes.size())));

  if (axes[0].n != (*vals).shape()[2])
    throw(SEPException(std::string("Axis 1 not the same (") +
                       std::to_string(axes[0].n) + std::string(",") +
                       std::to_string((*vals).shape()[2]) + std::string(")")));
  if (axes[1].n != (*vals).shape()[1])
    throw(SEPException(std::string("Axis 2 not the same (") +
                       std::to_string(axes[1].n) + std::string(",") +
                       std::to_string((*vals).shape()[1]) + std::string(")")));
  if (axes[2].n != (*vals).shape()[0])
    throw(SEPException(std::string("Axis 3 not the same (") +
                       std::to_string(axes[2].n) + std::string(",") +
                       std::to_string((*vals).shape()[0]) + std::string(")")));
  _mat.reset(new complex3D(boost::extents[axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (long long k = 0; k < axes[2].n; k++) {
    for (long long j = 0; j < axes[1].n; j++) {
      for (long long i = 0; i < axes[0].n; i++) {
        (*_mat)[k][j][i] = (*vals)[k][j][i];
      }
    }
  }
}
std::shared_ptr<complex3DReg> complex3DReg::window(
    const std::vector<int> &nw, const std::vector<int> &fw,
    const std::vector<int> &jw) const {
  const std::vector<SEP::axis> axes = getHyper()->getAxes();
  if (nw.size() != axes.size()) throw(SEPException("nw must of length 3"));
  if (fw.size() != axes.size()) throw(SEPException("fw must of length 3"));
  if (jw.size() != axes.size()) throw(SEPException("jw must of length 3"));
  std::vector<axis> aout;
  for (int i = 0; i < axes.size(); i++) {
    checkWindow(axes[i].n, nw[i], fw[i], jw[i]);
    aout.push_back(
        axis(nw[i], axes[i].o + axes[i].d * fw[i], axes[i].d * jw[i]));
  }
  std::shared_ptr<hypercube> hypOut(new hypercube(aout));
  std::shared_ptr<complex3DReg> out(new complex3DReg(hypOut));
  for (int i2 = 0; i2 < nw[2]; i2++) {
    for (int i1 = 0; i1 < nw[1]; i1++) {
      for (int i0 = 0; i0 < nw[0]; i0++) {
        (*out->_mat)[i2][i1][i0] =
            (*_mat)[fw[2] + i2 * jw[2]][fw[1] + i1 * jw[1]][fw[0] + i0 * jw[0]];
      }
    }
  }

  return out;
}