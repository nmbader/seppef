#include <byte1DReg.h>
using namespace SEP;
std::shared_ptr<byte1DReg> byte1DReg::clone() const {
  if (getSpaceOnly()) {
    std::shared_ptr<byte1DReg> x(new byte1DReg(getHyper()));

    return x;
  } else {
    std::shared_ptr<byte1DReg> x(new byte1DReg(getHyper(), *_mat));
    return x;
  }
}
std::shared_ptr<byte1DReg> byte1DReg::cloneSpace() const {
  std::shared_ptr<byte1DReg> x(new byte1DReg(getHyper()));
  x->_mat = 0;
  x->setSpace();

  return x;
}
byte1DReg::byte1DReg(const std::shared_ptr<byte6DReg> old, const int iax1,
                     const bool rev, const std::vector<int> &ipos,
                     const std::vector<int> &beg, const std::vector<int> &end) {
  std::vector<int> j(6, 1);
  std::vector<int> f(6, 0);
  std::vector<int> n(6, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }

  std::shared_ptr<byte6DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1));
  initNoData(hyperOut);
  int f1 = 0, j1 = 1;
  if (rev) {
    f1 = n[iax1] - 1;
    j1 = -1;
  }
  unsigned char *outv = getVals();
  unsigned char *inv = tmp->getVals();

  int ip1 = f1, i = 0;
  for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) outv[i] = inv[ip1];
}
byte1DReg::byte1DReg(const std::shared_ptr<byte5DReg> old, const int iax1,
                     const bool rev, const std::vector<int> &ipos,
                     const std::vector<int> &beg, const std::vector<int> &end) {
  std::vector<int> j(5, 1);
  std::vector<int> f(5, 0);
  std::vector<int> n(5, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }

  std::shared_ptr<byte5DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1));
  initNoData(hyperOut);
  int f1 = 0, j1 = 1;
  if (rev) {
    f1 = n[iax1] - 1;
    j1 = -1;
  }
  unsigned char *outv = getVals();
  unsigned char *inv = tmp->getVals();

  int ip1 = f1, i = 0;
  for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) outv[i] = inv[ip1];
}
byte1DReg::byte1DReg(const std::shared_ptr<byte4DReg> old, const int iax1,
                     const bool rev, const std::vector<int> &ipos,
                     const std::vector<int> &beg, const std::vector<int> &end) {
  std::vector<int> j(4, 1);
  std::vector<int> f(4, 0);
  std::vector<int> n(4, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }

  std::shared_ptr<byte4DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1));
  initNoData(hyperOut);
  int f1 = 0, j1 = 1;
  if (rev) {
    f1 = n[iax1] - 1;
    j1 = -1;
  }
  unsigned char *outv = getVals();
  unsigned char *inv = tmp->getVals();

  int ip1 = f1, i = 0;
  for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) outv[i] = inv[ip1];
}
byte1DReg::byte1DReg(const std::shared_ptr<byte3DReg> old, const int iax1,
                     const bool rev, const std::vector<int> &ipos,
                     const std::vector<int> &beg, const std::vector<int> &end) {
  std::vector<int> j(3, 1);
  std::vector<int> f(3, 0);
  std::vector<int> n(3, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }

  std::shared_ptr<byte3DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1));
  initNoData(hyperOut);
  int f1 = 0, j1 = 1;
  if (rev) {
    f1 = n[iax1] - 1;
    j1 = -1;
  }
  unsigned char *outv = getVals();
  unsigned char *inv = tmp->getVals();

  int ip1 = f1, i = 0;
  for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) outv[i] = inv[ip1];
}

byte1DReg::byte1DReg(const std::shared_ptr<byte2DReg> old, const int iax1,
                     const bool rev, const std::vector<int> &ipos,
                     const std::vector<int> &beg, const std::vector<int> &end) {
  std::vector<int> j(2, 1);
  std::vector<int> f(2, 0);
  std::vector<int> n(2, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }

  std::shared_ptr<byte2DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1));
  initNoData(hyperOut);
  int f1 = 0, j1 = 1;
  if (rev) {
    f1 = n[iax1] - 1;
    j1 = -1;
  }
  unsigned char *outv = getVals();
  unsigned char *inv = tmp->getVals();

  int ip1 = f1, i = 0;
  for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) outv[i] = inv[ip1];
}
byte1DReg::byte1DReg(const std::shared_ptr<byte1DReg> old, const int iax1,
                     const bool rev, const std::vector<int> &ipos,
                     const std::vector<int> &beg, const std::vector<int> &end) {
  std::vector<int> j(1, 1);
  std::vector<int> f(1, 0);
  std::vector<int> n(1, 1);
  std::vector<int> nd = old->getHyper()->getNs();
  // Figure out window
  for (auto i = 0; i < n.size(); i++) {
    f[i] = beg[i];
    if (iax1 == i)
      n[i] = end[i] - beg[i];
    else
      f[i] = ipos[i];
  }

  std::shared_ptr<byte1DReg> tmp = old->window(n, f, j);
  axis a1(n[iax1]);
  std::shared_ptr<hypercube> hyperOut(new hypercube(a1));
  initNoData(hyperOut);
  int f1 = 0, j1 = 1;
  if (rev) {
    f1 = n[iax1] - 1;
    j1 = -1;
  }
  unsigned char *outv = getVals();
  unsigned char *inv = tmp->getVals();

  int ip1 = f1, i = 0;
  for (auto i1 = 0; i1 < n[iax1]; i1++, i++, ip1 += j1) outv[i] = inv[ip1];
}

void byte1DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  setHyper(hyp);

  _vecType = "vec 1d int";
  const std::vector<SEP::axis> axes = hyp->getAxes();
  if (1 != axes.size()) throw(SEPException("axes must of length 2"));

  _mat.reset(new byte1D(boost::extents[axes[0].n]));
  setData(_mat->data());
}
void byte1DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                         const byte1D &vals) {
  setHyper(hyp);

  const std::vector<SEP::axis> axes = hyp->getAxes();

  if (axes.size() != 1)
    throw(SEPException(std::string("Axes size must be 1 is ") +
                       std::to_string(axes.size())));
  if (axes[0].n != vals.shape()[0])
    throw(SEPException(std::string("Axis 1 not the same (") +
                       std::to_string(axes[0].n) + std::string(",") +
                       std::to_string(vals.shape()[0]) + std::string(")")));
  _mat.reset(new byte1D(boost::extents[axes[0].n]));
  setData(_mat->data());
  for (long long i = 0; i < axes[0].n; i++) (*_mat)[i] = vals[i];
}
std::shared_ptr<byte1DReg> byte1DReg::window(const std::vector<int> &nw,
                                             const std::vector<int> &fw,
                                             const std::vector<int> &jw) const {
  const std::vector<SEP::axis> axes = getHyper()->getAxes();
  if (nw.size() != axes.size()) throw(SEPException("nw must of length 1"));
  if (fw.size() != axes.size()) throw(SEPException("fw must of length 1"));
  if (jw.size() != axes.size()) throw(SEPException("jw must of length 1"));
  std::vector<axis> aout;
  for (int i = 0; i < axes.size(); i++) {
    checkWindow(axes[i].n, nw[i], fw[i], jw[i]);
    aout.push_back(
        axis(nw[i], axes[i].o + axes[i].d * fw[i], axes[i].d * jw[i]));
  }
  std::shared_ptr<hypercube> hypOut(new hypercube(aout));
  std::shared_ptr<byte1DReg> out(new byte1DReg(hypOut));
  for (int i0 = 0; i0 < nw[0]; i0++) {
    (*out->_mat)[i0] = (*_mat)[fw[0] + i0 * jw[0]];
  }

  return out;
}
