#include <double2DReg.h>
using namespace SEP;

std::shared_ptr<double2DReg> double2DReg::clone() const {
  if (getSpaceOnly()) {
    std::shared_ptr<double2DReg> x(new double2DReg(getHyper()));

    return x;
  } else {
    std::shared_ptr<double2DReg> x(new double2DReg(getHyper(), *_mat));

    return x;
  }
}
std::shared_ptr<double2DReg> double2DReg::cloneSpace() const {
  std::shared_ptr<double2DReg> x(new double2DReg(getHyper()));
  x->_mat = 0;
  x->setSpace();

  return x;
}
double2DReg::double2DReg(const std::shared_ptr<double6DReg> old, const int iax1,
                         const bool rev1, const int iax2, const bool rev2,
                         const std::vector<int> &ipos,
                         const std::vector<int> &beg,
                         const std::vector<int> &end) {
  std::vector<int> j(6, 1);
  std::vector<int> f(6, 0);
  std::vector<int> n(6, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i || iax2 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }
  std::shared_ptr<double6DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]), a2(n[iax2]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1, a2));
  initNoData(hyperOut);
  int f1, j1, f2, j2;
  calcTraverse(n, iax1, rev1, f1, j1, iax2, rev2, f2, j2);
  int ip2 = f2, ip1 = f1, i = 0;
  double *outv = getVals();
  double *inv = tmp->getVals();
  for (auto i2 = 0; i2 < n[iax2]; i2++, ip2 += j2) {
    ip1 = f1;
    for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) {
      outv[i] = inv[ip1 + ip2];
    }
  }
}
double2DReg::double2DReg(const std::shared_ptr<double5DReg> old, const int iax1,
                         const bool rev1, const int iax2, const bool rev2,
                         const std::vector<int> &ipos,
                         const std::vector<int> &beg,
                         const std::vector<int> &end) {
  std::vector<int> j(5, 1);
  std::vector<int> f(5, 0);
  std::vector<int> n(5, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i || iax2 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }
  std::shared_ptr<double5DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]), a2(n[iax2]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1, a2));
  initNoData(hyperOut);
  int f1, j1, f2, j2;
  calcTraverse(n, iax1, rev1, f1, j1, iax2, rev2, f2, j2);
  int ip2 = f2, ip1 = f1, i = 0;
  double *outv = getVals();
  double *inv = tmp->getVals();
  for (auto i2 = 0; i2 < n[iax2]; i2++, ip2 += j2) {
    ip1 = f1;
    for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) {
      outv[i] = inv[ip1 + ip2];
    }
  }
}
double2DReg::double2DReg(const std::shared_ptr<double4DReg> old, const int iax1,
                         const bool rev1, const int iax2, const bool rev2,
                         const std::vector<int> &ipos,
                         const std::vector<int> &beg,
                         const std::vector<int> &end) {
  std::vector<int> j(4, 1);
  std::vector<int> f(4, 0);
  std::vector<int> n(4, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i || iax2 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }
  std::shared_ptr<double4DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]), a2(n[iax2]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1, a2));
  initNoData(hyperOut);
  int f1, j1, f2, j2;
  calcTraverse(n, iax1, rev1, f1, j1, iax2, rev2, f2, j2);
  int ip2 = f2, ip1 = f1, i = 0;
  double *outv = getVals();
  double *inv = tmp->getVals();
  for (auto i2 = 0; i2 < n[iax2]; i2++, ip2 += j2) {
    ip1 = f1;
    for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) {
      outv[i] = inv[ip1 + ip2];
    }
  }
}
double2DReg::double2DReg(const std::shared_ptr<double3DReg> old, const int iax1,
                         const bool rev1, const int iax2, const bool rev2,
                         const std::vector<int> &ipos,
                         const std::vector<int> &beg,
                         const std::vector<int> &end) {
  std::vector<int> j(3, 1);
  std::vector<int> f(3, 0);
  std::vector<int> n(3, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i || iax2 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }
  std::shared_ptr<double3DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]), a2(n[iax2]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1, a2));
  initNoData(hyperOut);
  int f1, j1, f2, j2;
  calcTraverse(n, iax1, rev1, f1, j1, iax2, rev2, f2, j2);
  int ip2 = f2, ip1 = f1, i = 0;
  double *outv = getVals();
  double *inv = tmp->getVals();
  for (auto i2 = 0; i2 < n[iax2]; i2++, ip2 += j2) {
    ip1 = f1;
    for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) {
      outv[i] = inv[ip1 + ip2];
    }
  }
}
double2DReg::double2DReg(const std::shared_ptr<double2DReg> old, const int iax1,
                         const bool rev1, const int iax2, const bool rev2,
                         const std::vector<int> &ipos,
                         const std::vector<int> &beg,
                         const std::vector<int> &end) {
  std::vector<int> j(2, 1);
  std::vector<int> f(2, 0);
  std::vector<int> n(2, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i || iax2 == i) n[i] = end[i] - beg[i];
  }
  std::shared_ptr<double2DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]), a2(n[iax2]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1, a2));
  initNoData(hyperOut);
  int f1, j1, f2, j2;
  calcTraverse(n, iax1, rev1, f1, j1, iax2, rev2, f2, j2);
  int ip2 = f2, ip1 = f1, i = 0;
  double *outv = getVals();
  double *inv = tmp->getVals();
  for (auto i2 = 0; i2 < n[iax2]; i2++, ip2 += j2) {
    ip1 = f1;
    for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) {
      outv[i] = inv[ip1 + ip2];
    }
  }
}
void double2DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);
  if (2 != axes.size()) throw(SEPException("must be 2-D hypercube"));

  _mat.reset(new double2D(boost::extents[axes[1].n][axes[0].n]));
  setData(_mat->data());
}
void double2DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                           const double2D &vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  _mat.reset(new double2D(boost::extents[axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (long long j = 0; j < axes[1].n; j++) {
    for (long long i = 0; i < axes[0].n; i++) {
      (*_mat)[j][i] = vals[j][i];
    }
  }
}
std::shared_ptr<double2DReg> double2DReg::window(
    const std::vector<int> &nw, const std::vector<int> &fw,
    const std::vector<int> &jw) const {
  const std::vector<SEP::axis> axes = getHyper()->getAxes();
  if (nw.size() != axes.size()) throw(SEPException("nw must of length 2"));
  if (fw.size() != axes.size()) throw(SEPException("fw must of length 2"));
  if (jw.size() != axes.size()) throw(SEPException("jw must of length 2"));
  std::vector<axis> aout;
  for (int i = 0; i < axes.size(); i++) {
    checkWindow(axes[i].n, nw[i], fw[i], jw[i]);
    aout.push_back(
        axis(nw[i], axes[i].o + axes[i].d * fw[i], axes[i].d * jw[i]));
  }
  std::shared_ptr<hypercube> hypOut(new hypercube(aout));
  std::shared_ptr<double2DReg> out(new double2DReg(hypOut));
  for (int i1 = 0; i1 < nw[1]; i1++) {
    for (int i0 = 0; i0 < nw[0]; i0++) {
      (*out->_mat)[i1][i0] = (*_mat)[fw[1] + i1 * jw[1]][fw[0] + i0 * jw[0]];
    }
  }

  return out;
}