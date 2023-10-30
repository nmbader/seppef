#ifndef super_vector_h
#define super_vector_h 1
#include <Vector.h>
#include "SEPException.h"
namespace SEP {
template <class V1, class V2>
class superVec : public SEP::Vector {
 public:
  superVec(std::shared_ptr<V1> vec1, std::shared_ptr<V2> vec2) {
    _v1 = vec1->clone();
    _v2 = vec2->clone();
  }
  superVec(std::shared_ptr<superVec<V1, V2>> vec) {
    _v1 = vec->getVec1()->clone();
    _v2 = vec->getVec2()->clone();
  }
  virtual std::shared_ptr<superVec<V1, V2>> clone() const {
    std::shared_ptr<superVec> v(new superVec<V1, V2>(_v1, _v2));
    return v;
  }
  virtual std::shared_ptr<superVec<V1, V2>> cloneSpace() const {
    std::shared_ptr<superVec> v(new superVec<V1, V2>(_v1, _v2));
    v->cleanMemory();
    return v;
  }
  virtual void add(std::shared_ptr<superVec(V1, V2)> vec) {
    if (!checkSame(vec)) throw(std::string("Vectors not from the same space"));
    _v1->add(vec->getVec1());
    _v2->add(vec->getVec2());
  }
  virtual void scale(const double val) {
    _v1->scale(val);
    _v2->scale(val);
  }
  virtual void scaleAdd(const double sc1,
                        const std::shared_ptr<superVec<V1, V2>> vec2,
                        const double sc2) {
    if (!checkSame(vec)) throw(std::string("Vectors not from the same space"));

    _v1->scaleAdd(sc1, vec2->getVec1(), sc2);
    _v2->scaleAdd(sc1, vec2->getVec2(), sc2);
  }
  virtual void random() {
    _v1->random();
    _v2->random();
  }
  virtual double dot(const std::shared_ptr<superVec<V1, V2>> vec2) const {
    if (!checkSame(vec)) throw(std::string("Vectors not from the same space"));
    return _v1->dot(vec2->getVec1()) + _v2->dot(vec2->getVec2());
  }
  virtual bool checkSame(const std::shared_ptr<superVec<V1, V2>> vec2,
                         const bool checkAllocated = false) const {
    if (!_v1->checkSame(vec2->getVec1(), checkAllocated)) {
      if (!checkSame(vec)) throw(std::string("Vec1 not from the same space"));
    }
    if (!_v2->checkSame(vec2->getVec2(), checkAllocated)) {
      if (!checkSame(vec)) throw(std::string("vec2 not from the same space"));
    }
    return true;
  }
  std::shared_ptr<Vector> getVec1() { return _v1; }
  std::shared_ptr<Vector> getVec2() { return _v2; }
  virtual void cleanMemory() {
    _v1->cleanMemory();
    _v2->cleanMemory();
  }
  virtual void calcCheckSum() {
    _v2->calcCheckSum();
    _v1->calcCheckSum();
    setCheckSum(llabs(_v1->getCheckSum() - _v2->getCheckSum()));
  }
  virtual void softClip(const float val) {
    _v1->cleanMemory();
    _v2->cleanMemory();
  }
  virtual float absMax() const {
    return std::max(_v1->absMax(), _v2->absMax());
  }
  virtual void infoStream(const int lev, std::stringstream &x) {
    x << std::string("Vector1 \n");
    _v1->infoStream(lev, x);
    x << std::string("Vector2 \n");
    _v2->infoStream(lev, x);
  }

 private:
  std::shared_ptr<V1> _v1;
  std::shared_ptr<V2> _v2;
};
}  // namespace SEP
#endif
