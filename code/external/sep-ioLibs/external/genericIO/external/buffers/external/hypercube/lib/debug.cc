#include "debug.h"
#include "SEPException.h"
using namespace SEP;
void debug::addLine(std::stringstream &x) {
  _lines.push_back(x.str());
  if (_lines.size() >= _nflush) flush();
}

void errDebug::flush() {
  for (int i = 0; i < _lines.size(); i++) std::cerr << _lines[i] << std::endl;
  _lines.resize(0);
}
void fileDebug::flush() {
  std::string x;
  for (int i = 0; i < _lines.size(); i++) x = x + _lines[i] + "\n";
  _file << x;
  _file.flush();
  _lines.resize(0);
}
void fileDebug::close() { _file.close(); }
fileDebug::fileDebug(std::string &name, int nflush) {
  _file.open(name);
  if (!_file) SEPException("Trouble opening file");
  setFlush(nflush);
}

std::shared_ptr<debug> debugging::getDebug(const std::string &name) {
  if (_debugs.count(name) == 0)
    throw SEPException(std::string("Request for unknown debug"));
  return _debugs[name];
}

std::vector<std::string> getDebug(const std::string &buffer) {
  debugging *mydebug = &debugging::instance();
  return mydebug->getDebug(buffer)->returnLines();
}
