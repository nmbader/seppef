#include <float6DReg.h>

using namespace SEP;

std::shared_ptr<float6DReg> float6DReg::clone() const {
  if (getSpaceOnly()) {
    std::shared_ptr<float6DReg> x(new float6DReg(getHyper()));

    return x;
  } else {
    std::shared_ptr<float6DReg> x(new float6DReg(getHyper(), *_mat));

    return x;
  }
}
std::shared_ptr<float6DReg> float6DReg::cloneSpace() const {
  std::shared_ptr<float6DReg> x(new float6DReg(getHyper()));
  x->_mat = 0;

  x->setSpace();
  return x;
}

void float6DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  if (axes.size() != 6)
    throw(SEPException(std::string("Axes size must be 6 is ") +
                       std::to_string(axes.size())));

  _mat.reset(new float6D(boost::extents[axes[5].n][axes[4].n][axes[3].n]
                                       [axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
}
void float6DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                          const float6D &vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  if (axes.size() != 6)
    throw(SEPException(std::string("Axes size must be 6 is ") +
                       std::to_string(axes.size())));

  if (axes[0].n != vals.shape()[5])
    throw(SEPException(std::string("Axis 1 not the same (") +
                       std::to_string(axes[0].n) + std::string(",") +
                       std::to_string(vals.shape()[5]) + std::string(")")));
  if (axes[1].n != vals.shape()[4])
    throw(SEPException(std::string("Axis 2 not the same (") +
                       std::to_string(axes[1].n) + std::string(",") +
                       std::to_string(vals.shape()[4]) + std::string(")")));
  if (axes[2].n != vals.shape()[3])
    throw(SEPException(std::string("Axis 3 not the same (") +
                       std::to_string(axes[2].n) + std::string(",") +
                       std::to_string(vals.shape()[3]) + std::string(")")));
  if (axes[3].n != vals.shape()[2])
    throw(SEPException(std::string("Axis 4 not the same (") +
                       std::to_string(axes[3].n) + std::string(",") +
                       std::to_string(vals.shape()[2]) + std::string(")")));
  if (axes[4].n != vals.shape()[1])
    throw(SEPException(std::string("Axis 5 not the same (") +
                       std::to_string(axes[4].n) + std::string(",") +
                       std::to_string(vals.shape()[1]) + std::string(")")));
  if (axes[5].n != vals.shape()[0])
    throw(SEPException(std::string("Axis 6 not the same (") +
                       std::to_string(axes[5].n) + std::string(",") +
                       std::to_string(vals.shape()[0]) + std::string(")")));
  _mat.reset(new float6D(boost::extents[axes[5].n][axes[4].n][axes[3].n]
                                       [axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (size_t n = 0; n < axes[5].n; n++) {
    for (size_t m = 0; m < axes[4].n; m++) {
      for (long long l = 0; l < axes[3].n; l++) {
        for (long long k = 0; k < axes[2].n; k++) {
          for (long long j = 0; j < axes[1].n; j++) {
            for (long long i = 0; i < axes[0].n; i++) {
              (*_mat)[n][m][l][k][j][i] = vals[n][m][l][k][j][i];
            }
          }
        }
      }
    }
  }
}

std::shared_ptr<float6DReg> float6DReg::window(
    const std::vector<int> &nw, const std::vector<int> &fw,
    const std::vector<int> &jw) const {
  const std::vector<SEP::axis> axes = getHyper()->getAxes();
  if (nw.size() != axes.size()) throw(SEPException("nw must of length 6"));
  if (fw.size() != axes.size()) throw(SEPException("fw must of length 6"));
  if (jw.size() != axes.size()) throw(SEPException("jw must of length 6"));
  std::vector<axis> aout;
  for (int i = 0; i < axes.size(); i++) {
    checkWindow(axes[i].n, nw[i], fw[i], jw[i]);
    aout.push_back(
        axis(nw[i], axes[i].o + axes[i].d * fw[i], axes[i].d * jw[i]));
  }
  std::shared_ptr<hypercube> hypOut(new hypercube(aout));
  std::shared_ptr<float6DReg> out(new float6DReg(hypOut));
  for (int i5 = 0; i5 < nw[5]; i5++) {
    for (int i4 = 0; i4 < nw[4]; i4++) {
      for (int i3 = 0; i3 < nw[3]; i3++) {
        for (int i2 = 0; i2 < nw[2]; i2++) {
          for (int i1 = 0; i1 < nw[1]; i1++) {
            for (int i0 = 0; i0 < nw[0]; i0++) {
              (*out->_mat)[i5][i4][i3][i2][i1][i0] =
                  (*_mat)[fw[5] + i5 * jw[5]][fw[4] + i4 * jw[4]]
                         [fw[3] + i3 * jw[3]][fw[2] + i2 * jw[2]]
                         [fw[1] + i1 * jw[1]][fw[0] + i0 * jw[0]];
            }
          }
        }
      }
    }
  }
  return out;
}
