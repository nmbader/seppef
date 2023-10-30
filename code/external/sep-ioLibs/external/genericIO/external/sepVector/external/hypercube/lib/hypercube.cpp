#include <SEPException.h>
#include <hypercube.h>
#include "math.h"
using namespace SEP;
hypercube::hypercube(const hypercube &hyper) {
  std::vector<axis> axes;
  axes = hyper.getAxes();
  initNd(axes);
}
hypercube::hypercube(const std::shared_ptr<hypercube> hyper) {
  std::vector<axis> axes;
  axes = hyper->getAxes();
  initNd(axes);
}
hypercube::hypercube(const int ndim) {
  std::vector<axis> axes;
  for (int i = 0; i < ndim; i++) axes.push_back(axis(1));
  initNd(axes);
}
hypercube::hypercube(const std::vector<axis> &axes) { initNd(axes); }
void hypercube::initNd(const std::vector<axis> &ax) {
  this->axes.resize(0);
  this->n123 = 1;
  for (int i = 0; i < ax.size(); i++) {
    axes.push_back(ax[i]);
    this->n123 = this->n123 * (long long)ax[i].n;
  }
}
std::vector<axis> hypercube::getAxes(const int nmin) const {
  std::vector<axis> ax = getAxes();
  for (int i = ax.size(); i < nmin; i++) ax.push_back(axis(1));
  return ax;
}
bool hypercube::checkSame(const std::shared_ptr<hypercube> hyper2) const {
  if (hyper2->getAxes().size() != axes.size()){
    return false;
    throw SEPException("Axes not the same length");
  }
  for (int i = 0; i < axes.size(); i++) {
    const axis a = hyper2->getAxis(i + 1);
    if (a.n != axes[i].n){
    return false;

      throw SEPException(std::string("Axis ") + std::to_string(i + 1) +
                         std::string(" n do not match ") +
                         std::to_string(axes[i].n) + std::string(" ") +
                         std::to_string(a.n));
    }
    if (fabs((a.d - axes[i].d) / a.d) > 1e-3){
      return false;
      throw SEPException(std::string("Axis ") + std::to_string(i + 1) +
                         std::string(" d do not match ") +
                         std::to_string(axes[i].d) + std::string(" ") +
                         std::to_string(a.d));
    }
    if (fabs((a.o - axes[i].o) / a.d) > 1e-3){
      return false;
      throw SEPException(std::string("Axis ") + std::to_string(i + 1) +
                         std::string(" o do not match  ") +
                         std::to_string(axes[i].o) + std::string(" ") +
                         std::to_string(a.o));
    }
  }
  return true;
}

void hypercube::setAxis(const int idim, const axis &myaxis) {
  if (idim < 1 || idim > axes.size()) {
    throw SEPException(std::string("idim=") + std::to_string(idim) +
                       std::string(" axes.size()=") +
                       std::to_string(axes.size()));
  }
  this->axes[idim - 1] = myaxis;
}
void hypercube::setAxes(const std::vector<axis> &axes) { this->initNd(axes); }
std::vector<axis> hypercube::returnAxes(const int nmax) const {
  std::vector<axis> ax;
  for (int i = 0; i < nmax; i++) {
    if (i + 1 <= this->axes.size()) {
      ax.push_back(axes[i]);
    } else {
      ax.push_back(axis(1));
    }
  }
  return ax;
}
void hypercube::infoStream(std::stringstream &x) {
  for (int i = 0; i < (int)axes.size(); i++) {
    x << "Axis " << std::to_string(i + 1);
    axes[i].infoStream(x);
    x << "\n";
  }
}
axis hypercube::getAxis(const int idim) const {
  if (idim < 1 || idim > this->axes.size()) {
    throw SEPException(std::string("IDIM=") + std::to_string(idim) +
                       std::string(" axes.size()=") +
                       std::to_string(axes.size()));
  }
  axis myaxis = this->axes[idim - 1];
  return myaxis;
}
int hypercube::getNdimG1() const {
  int nd = this->axes.size();
  for (int i = this->axes.size() - 1; i >= 0; i--) {
    if (axes[i].n > 1) return nd;
    nd--;
  }

  return nd;
}

std::vector<axis> hypercube::getAxes() const { return returnAxes(axes.size()); }
std::vector<int> hypercube::getNs() const {
  int i;
  std::vector<int> n;
  for (i = 0; i < this->axes.size(); i++) n.push_back(axes[i].n);
  return n;
}
bool hypercube::sameSize(const hypercube &other) const {
  if (this->getNdim() != other.getNdim()) return false;
  for (int i = 0; i < this->getNdim(); i++) {
    axis ax1 = this->getAxis(i + 1);
    axis ax2 = other.getAxis(i + 1);
    if (ax1.n != ax2.n) return false;
  }
  return true;
}
bool hypercube::sameSize(const std::shared_ptr<hypercube> &other) const {
  if (this->getNdim() != other->getNdim()) return false;
  for (int i = 0; i < this->getNdim(); i++) {
    axis ax1 = this->getAxis(i + 1);
    axis ax2 = other->getAxis(i + 1);
    if (ax1.n != ax2.n) return false;
  }
  return true;
}
