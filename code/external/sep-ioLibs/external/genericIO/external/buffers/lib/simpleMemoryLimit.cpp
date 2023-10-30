#include "simpleMemoryLimit.h"
#include <iostream>
using namespace SEP::IO;

simpleMemoryLimit::simpleMemoryLimit(const size_t cleanAt) {
  _maxMem = cleanAt;
}
void simpleMemoryLimit::updateRecentBuffers(const std::vector<int> &bufs) {
  for (auto i = 0; i < bufs.size(); i++) {
    _recent[bufs[i]] = _ibuf;
    _status[bufs[i]] = 0;
  }
  _ibuf += 1;
}
std::shared_ptr<memoryReduce> simpleMemoryLimit::changeBufferState(
    const long memChange) {
  _curMem += memChange;
  //  std::cerr << "Current usage " << _curMem / 1024 / 1024 << " "
  //        << memChange / 1024 / 1024 << " allowed=" << _maxMem / 1024 / 1024
  //    << std::endl;
  std::vector<int> comp, disk;
  if (_curMem > _maxMem) {
    int iold = _compressed;
    _compressed = std::min(_ibuf - 1, (size_t)(_compressed + 2));
    for (auto i = _recent.begin(); i != _recent.end(); ++i) {
      if (i->second <= _compressed && (int)i->second > iold) {
        if (_status[i->first] == 0) {
          comp.push_back(i->first);
          _status[i->first] = 1;
        }
      }
    }

    if (comp.size() == 0) {
      int iold = _toDisk;
      _toDisk = std::min(_ibuf - 1, (size_t)_toDisk + 2);
      for (auto i = _recent.begin(); i != _recent.end(); ++i) {
        if (i->second <= _toDisk && (int)i->second > iold) {
          if (_status[i->first] != 2) {
            disk.push_back(i->first);
            _status[i->first] = 2;
          }
        }
      }
    }
  }
  std::shared_ptr<memoryReduce> x(new memoryReduce(comp, disk));
  return x;
}
