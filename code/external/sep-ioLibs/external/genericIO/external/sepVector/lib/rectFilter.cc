#include "rectFilter.h"
using namespace SEP;

rectFilter1D::rectFilter1D(const std::vector<int>& box,
                           const std::vector<int>& f0, bool pef) {
  setBasics(box, f0, pef);
  axis a(box[0]);
  std::shared_ptr<hypercube> hyper(new hypercube(a));

  initNoData(hyper);

  zeroNonCoefs();
}
std::shared_ptr<rectFilter1D> rectFilter1D::clone() const {
  std::shared_ptr<rectFilter1D> f(new rectFilter1D(_n, _f, _pef));
  for (auto i1 = 0; i1 < getHyper()->getAxis(1).n; i1++) {
    (*f->_mat)[i1] = (*_mat)[i1];
  }
  return f;
}
std::shared_ptr<rectFilter1D> rectFilter1D::cloneSpace() const {
  std::shared_ptr<rectFilter1D> f(new rectFilter1D(_n, _f, _pef));
  f->setSpace();
  return f;
}

void rectFilter1D::zeroNonCoefs() {
  for (auto i1 = 0; i1 < _f[0]; i1++) {
    (*_mat)[i1] = 0;
  }
}

void rectFilter2D::zeroNonCoefs() {
  for (auto i2 = 0; i2 <= _f[1]; i2++) {
    int end = _n[0];
    if (_f[1] == i2) end = _f[0] + 1;
    for (auto i1 = 0; i1 < end; i1++) {
      (*_mat)[i2][i1] = 0;
    }
  }
}

rectFilter2D::rectFilter2D(const std::vector<int>& box,
                           const std::vector<int>& f0, bool pef) {
  setBasics(box, f0, pef);
  std::shared_ptr<hypercube> hyper(new hypercube(box[0], box[1]));

  initNoData(hyper);

  zero();
}
std::shared_ptr<rectFilter2D> rectFilter2D::clone() const {
  std::shared_ptr<rectFilter2D> f(new rectFilter2D(_n, _f, _pef));

  for (auto i2 = 0; i2 < getHyper()->getAxis(2).n; i2++) {
    for (auto i1 = 0; i1 < getHyper()->getAxis(1).n; i1++) {
      (*f->_mat)[i2][i1] = (*_mat)[i2][i1];
    }
  }
  return f;
}
std::shared_ptr<rectFilter2D> rectFilter2D::cloneSpace() const {
  std::shared_ptr<rectFilter2D> f(new rectFilter2D(_n, _f, _pef));
  f->setSpace();
  return f;
}

rectFilter3D::rectFilter3D(const std::vector<int>& box,
                           const std::vector<int>& f0, bool pef) {
  setBasics(box, f0, pef);
  std::shared_ptr<hypercube> hyper(new hypercube(box[0], box[1], box[2]));
  initNoData(hyper);
}
std::shared_ptr<rectFilter3D> rectFilter3D::clone() const {
  std::shared_ptr<rectFilter3D> f(new rectFilter3D(_n, _f));
  for (auto i3 = 0; i3 < getHyper()->getAxis(3).n; i3++) {
    for (auto i2 = 0; i2 < getHyper()->getAxis(2).n; i2++) {
      for (auto i1 = 0; i1 < getHyper()->getAxis(1).n; i1++) {
        (*f->_mat)[i3][i2][i1] = (*_mat)[i3][i2][i1];
      }
    }
  }
  return f;
}
std::shared_ptr<rectFilter3D> rectFilter3D::cloneSpace() const {
  std::shared_ptr<hypercube> hyper(new hypercube(getHyper()->getN123()));
  std::shared_ptr<rectFilter3D> f(new rectFilter3D(_n, _f));
  f->setSpace();
  return f;
}
void rectFilter3D::zeroNonCoefs() {
  for (auto i3 = 0; i3 <= _f[2]; i3++) {
    for (auto i2 = 0; i2 <= _f[1]; i2++) {
      int end = _n[0];
      if (_f[1] == i2) end = _f[0] + 1;
      for (auto i1 = 0; i1 < end; i1++) {
        (*_mat)[i3][i2][i1] = 0;
      }
    }
  }
}

void rectFilter4D::zeroNonCoefs() {
  for (auto i4 = 0; i4 <= _f[3]; i4++) {
    for (auto i3 = 0; i3 <= _f[2]; i3++) {
      for (auto i2 = 0; i2 <= _f[1]; i2++) {
        int end = _n[0];
        if (_f[1] == i2) end = _f[0] + 1;
        for (auto i1 = 0; i1 < end; i1++) {
          (*_mat)[i4][i3][i2][i1] = 0;
        }
      }
    }
  }
}

rectFilter4D::rectFilter4D(const std::vector<int>& box,
                           const std::vector<int>& f0, bool pef) {
  setBasics(box, f0, pef);

  std::shared_ptr<hypercube> hyper(
      new hypercube(box[0], box[1], box[2], box[3]));
  initNoData(hyper);
  zero();
}
std::shared_ptr<rectFilter4D> rectFilter4D::clone() const {
  std::shared_ptr<rectFilter4D> f(new rectFilter4D(_n, _f));
  for (auto i4 = 0; i4 < getHyper()->getAxis(4).n; i4++) {
    for (auto i3 = 0; i3 < getHyper()->getAxis(3).n; i3++) {
      for (auto i2 = 0; i2 < getHyper()->getAxis(2).n; i2++) {
        for (auto i1 = 0; i1 < getHyper()->getAxis(1).n; i1++) {
          (*f->_mat)[i4][i3][i2][i1] = (*_mat)[i4][i3][i2][i1];
        }
      }
    }
  }
  return f;
}
std::shared_ptr<rectFilter4D> rectFilter4D::cloneSpace() const {
  std::shared_ptr<rectFilter4D> f(new rectFilter4D(_n, _f));

  f->setSpace();
  return f;
}
