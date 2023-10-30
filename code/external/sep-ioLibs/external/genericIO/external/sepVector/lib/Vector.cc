#include <Vector.h>
using namespace SEP;
std::string Vector::info(const std::string& nm, const int lev) {
  std::stringstream ss;
  ss << nm << std::endl;
  ss << "---------------------------" << std::endl;
  infoStream(lev, ss);
  ss << "---------------------------" << std::endl;
  return ss.str();
}
