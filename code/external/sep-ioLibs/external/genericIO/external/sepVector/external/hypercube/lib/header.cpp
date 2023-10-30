#include "header.h"
#include <cstring>
#include "SEPException.h"
using namespace SEP;

header::header(const std::vector<key> keys, const int nh) {
  _keys = keys;
  int nsz = 0;

  for (auto i = 0; i < keys.size(); i++) {
    _key_offset[keys[i].name()] = nsz;
    nsz += getDataTypeSize(keys[i].type());
    _key_index[keys[i].name()] = i;
  }
  std::vector<unsigned char> head(nsz);
  _head.resize(nh, head);
}

std::shared_ptr<header> header::clone() {
  std::shared_ptr<header> x(new header(_keys, _head.size()));
  x->setHeaders(_head);
  return x;
}
void header::addKey(const std::string name, const dataType typ) {
  _key_index[name] = _keys.size();
  _keys.push_back(key(name, typ));
  _key_offset[name] = _nsz;
  size_t nsz = _nsz + getDataTypeSize(typ);

  if (_head.size() == 0) {
    for (size_t ih = 0; ih < _head.size(); ih++) {
      _head[ih].resize(nsz);
    }
  }

  _nsz = nsz;
}
void header::setFloatKey(const std::string name,
                         const std::vector<float> vals) {
  size_t off = getOffset(name);
  if (vals.size() != _head.size())
    throw(SEPException(
        std::string("Size of vals not equal to number of headers")));
  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(_head[ih].data() + off, vals.data() + sizeof(float) * ih,
           sizeof(float));
}
void header::setIntKey(const std::string name, const std::vector<int> vals) {
  size_t off = getOffset(name);
  if (vals.size() != _head.size())
    throw(SEPException(
        std::string("Size of vals not equal to number of headers")));
  if (_keys[getIndex(name)].type() != DATA_INT)
    throw(SEPException(std::string("Key type is not integer")));

  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(_head[ih].data() + off, vals.data() + sizeof(int) * ih, sizeof(int));
}
void header::setDoubleKey(const std::string name,
                          const std::vector<double> vals) {
  size_t off = getOffset(name);
  if (vals.size() != _head.size())
    throw(SEPException(
        std::string("Size of vals not equal to number of headers")));
  if (_keys[getIndex(name)].type() != DATA_DOUBLE)
    throw(SEPException(std::string("Key type is not double")));
  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(_head[ih].data() + off, vals.data() + sizeof(double) * ih,
           sizeof(double));
}
void header::setShortKey(const std::string name,
                         const std::vector<short> vals) {
  size_t off = getOffset(name);
  if (vals.size() != _head.size())
    throw(SEPException(
        std::string("Size of vals not equal to number of headers")));
  if (_keys[getIndex(name)].type() != DATA_SHORT)
    throw(SEPException(std::string("Key type is not short")));
  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(_head[ih].data() + off, vals.data() + sizeof(short) * ih,
           sizeof(short));
}

void header::setHeaders(std::vector<std::vector<unsigned char>>& head) {
  if (head[0].size() != _head.size())
    throw(
        SEPException(std::string("Header size not equal to number of keys ")));
  if (_head.size() != 0 && head.size() != _head.size())
    throw SEPException(
        std::string("Number of headers not equal to number set"));

  _head = head;
}
std::vector<float> header::getFloatKey(const std::string name) const {
  size_t off = getOffset(name);
  std::vector<float> vals(_head.size());
  if (_keys[getIndex(name)].type() != DATA_FLOAT)
    throw(SEPException(std::string("Key type is not float")));
  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(&vals[ih], _head[ih].data() + off, sizeof(float));
  return vals;
}
std::vector<double> header::getDoubleKey(const std::string name) const {
  size_t off = getOffset(name);
  std::vector<double> vals(_head.size());
  if (_keys[getIndex(name)].type() != DATA_DOUBLE)
    throw(SEPException(std::string("Key type is not double")));

  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(&vals[ih], _head[ih].data() + off, sizeof(double));
  return vals;
}
std::vector<int> header::getIntKey(const std::string name) const {
  size_t off = getOffset(name);
  std::vector<int> vals(_head.size());
  if (_keys[getIndex(name)].type() != DATA_INT)
    throw(SEPException(std::string("Key type is not int")));

  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(&vals[ih], _head[ih].data() + off, sizeof(int));
  return vals;
}
std::vector<short> header::getShortKey(const std::string name) const {
  size_t off = getOffset(name);
  std::vector<short> vals(_head.size());
  if (_keys[getIndex(name)].type() != DATA_SHORT)
    throw(SEPException(std::string("Key type is not short")));

  for (size_t ih = 0; ih < _head.size(); ih++)
    memcpy(&vals[ih], _head[ih].data() + off, sizeof(short));
  return vals;
}
float header::getFloatKeyVal(const std::string name, const int index) const {
  size_t off = getOffset(name);
  float val;
  if (_keys[getIndex(name)].type() != DATA_FLOAT)
    throw(SEPException(std::string("Key type is not float")));

  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));

  memcpy(&val, _head[index].data() + off, sizeof(float));
  return val;
}
double header::getDoubleKeyVal(const std::string name, const int index) const {
  size_t off = getOffset(name);
  double val;
  if (_keys[getIndex(name)].type() != DATA_DOUBLE)
    throw(SEPException(std::string("Key type is not double")));
  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));

  memcpy(&val, _head[index].data() + off, sizeof(double));
  return val;
}
int header::getIntKeyVal(const std::string name, const int index) const {
  size_t off = getOffset(name);
  int val;
  if (_keys[getIndex(name)].type() != DATA_INT)
    throw(SEPException(std::string("Key type is not int")));
  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));

  memcpy(&val, _head[index].data() + off, sizeof(int));
  return val;
}
short header::getShortKeyVal(const std::string name, const int index) const {
  size_t off = getOffset(name);
  short val;
  if (_keys[getIndex(name)].type() != DATA_SHORT)
    throw(SEPException(std::string("Key type is not short")));
  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));

  memcpy(&val, _head[index].data() + off, sizeof(short));
  return val;
}
void header::setFloatKeyVal(const std::string name, const int index,
                            const float val) {
  size_t off = getOffset(name);
  if (_keys[getIndex(name)].type() != DATA_FLOAT)
    throw(SEPException(std::string("Key type is not float")));
  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));

  memcpy(_head[index].data() + off, &val, sizeof(float));
}
void header::setShortKeyVal(const std::string name, const int index,
                            const short val) {
  size_t off = getOffset(name);
  if (_keys[getIndex(name)].type() != DATA_SHORT)
    throw(SEPException(std::string("Key type is not short")));
  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));
  memcpy(_head[index].data() + off, &val, sizeof(short));
}
void header::setDoubleKeyVal(const std::string name, const int index,
                             const double val) {
  size_t off = getOffset(name);
  if (_keys[getIndex(name)].type() != DATA_DOUBLE)
    throw(SEPException(std::string("Key type is not double")));
  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));
  memcpy(_head[index].data() + off, &val, sizeof(double));
}

void header::setIntKeyVal(const std::string name, const int index,
                          const int val) {
  size_t off = getOffset(name);
  if (_keys[getIndex(name)].type() != DATA_INT)
    throw(SEPException(std::string("Key type is not int")));
  if (index < 0 || index >= _head.size())
    throw SEPException(
        std::string("Header requested out of range must be between 0 and " +
                    std::to_string(_head.size())));
  memcpy(_head[index].data() + off, &val, sizeof(int));
}
