#include <axis.h>
#include <math.h>
#include "SEPException.h"

using std::string;
using namespace SEP;

axis::axis(const int n, float o, float d, std::string label, std::string unit) {
  this->n = n;
  this->o = o;
  this->d = d;
  this->label = label;
  this->unit = unit;
}

bool axis::same_axis(const axis& ax) const {
  bool match = true;
  if (fabs(this->d) < 1e-30)
    throw(SEPException(std::string("d to small ") + std::to_string(this->d)));
  if (this->n != ax.n) match = false;
  if (fabs((this->o - ax.o) / this->d) > .01) match = false;
  if (fabs(this->d / ax.d - 1.) > .01) match = false;
  return match;
}

void axis::infoStream(std::stringstream& stream) {
  stream << std::string(" n=") << std::to_string(n) << std::string(" o=")
         << std::to_string(o);
  stream << std::string("d=") << std::to_string(d) << std::string(" label=")
         << label;
}
