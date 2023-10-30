#include <superVec.h>
#include <cmath>
#include <iostream>
using namespace SEP;
superVec::superVec(std::shared_ptr<Vector> vec1, std::shared_ptr<Vector> vec2) {
  _v1 = vec1->clone();
  _v2 = vec2->clone();
}
superVec::superVec(std::shared_ptr<superVec> s) {
  _v1 = s->getVec1()->clone();
  _v2 = s->getVec2()->clone();
}
std::shared_ptr<superVec<V1, V2>> superVec::clone() const {
  std::shared_ptr<superVec<V1, V2>> v(new superVec(_v1, _v2));
  return v;
}
std::shared_ptr<superVec<V1, V2>> superVec::cloneSpace() const {
  std::shared_ptr<superVec<V1, V2>> v(new superVec(_v1, _v2));
  v->cleanMemory();
  return v;
}
void superVec::calcCheckSum() {
  _v2->calcCheckSum();
  _v1->calcCheckSum();
  setCheckSum(llabs(_v1->getCheckSum() - _v2->getCheckSum()));
}
void superVec::add(const std::shared_ptr<Vector> vec) {
  if (!checkSame(vec)) throw(std::string("Vectors not from the same space"));
  std::shared_ptr<superVec> vec2H = std::dynamic_pointer_cast<superVec>(vec);
  _v1->add(vec2H->getVec1());
  _v2->add(vec2H->getVec2());
}
void superVec::scale(const double val) {
  _v1->scale(val);
  _v2->scale(val);
}
void superVec::scaleAdd(const double sc1, const std::shared_ptr<Vector> vec2,
                        const double sc2) {
  if (!checkSame(vec2)) throw(std::string("Vectors not from the same space"));
  std::shared_ptr<superVec> vec2H = std::dynamic_pointer_cast<superVec>(vec2);
  _v1->scaleAdd(sc1, vec2H->getVec1(), sc2);
  _v2->scaleAdd(sc1, vec2H->getVec2(), sc2);
}
void superVec::random() {
  _v1->random();
  _v2->random();
}
void superVec::softClip(const float val) {
  _v1->softClip(val);
  _v2->softClip(val);
}
float superVec::absMax() const {
  return std::max(_v1->absMax(), _v2->absMax());
}

double superVec::dot(const std::shared_ptr<Vector> vec2) const {
  if (!checkSame(vec2)) throw(std::string("Vectors not from the same space"));
  std::shared_ptr<superVec> vec2H = std::dynamic_pointer_cast<superVec>(vec2);
  return _v1->dot(vec2H->getVec1()) + _v2->dot(vec2H->getVec2());
}
void superVec::cleanMemory() {
  _v1->cleanMemory();
  _v2->cleanMemory();
}
void superVec::infoStream(const int lev, std::stringstream &x) {
  x << std::string("Vector1 \n");
  _v1->infoStream(lev, x);
  x << std::string("Vector2 \n");
  _v2->infoStream(lev, x);
}

bool superVec::checkSame(const std::shared_ptr<Vector> vec2,
                         const bool checkAllocated) const {
  std::shared_ptr<superVec> vec2H = std::dynamic_pointer_cast<superVec>(vec2);
  if (!vec2H) {
    std::cerr << "Not a superVec" << std::endl;
    return false;
  }
  if (!_v1->checkSame(vec2H->getVec1(), checkAllocated)) {
    std::cerr << "Vec1 does not match" << std::endl;
    return false;
  }
  if (!_v2->checkSame(vec2H->getVec2(), checkAllocated)) {
    std::cerr << "Vec2 does not match" << std::endl;
    return false;
  }
  return true;
}
