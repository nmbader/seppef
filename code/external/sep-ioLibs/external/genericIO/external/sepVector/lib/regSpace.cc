
#include <iostream>

#include "regSpace.h"
using namespace SEP;
void regSpace::checkWindow(const int n, const int nw, const int fw,
                           const int jw) const {
  if (fw < 0)
    throw(SEPException(std::string("fw must be 0 or greater=") +
                       std::to_string(fw)));
  if (jw < 1)
    throw(SEPException(std::string("jw must be 1 or greater=") +
                       std::to_string(jw)));
  if (nw < 1)
    throw(SEPException(std::string("nw must be 1 or greater=") +
                       std::to_string(nw)));
  if (fw + jw * (nw - 1) >= n)
    throw(SEPException(std::string("fw+jw*(nw-1) <n ") + std::to_string(fw) +
                       std::string("=fw jw=") + std::to_string(jw) +
                       std::string(" nw=") + std::to_string(nw) +
                       std::string(" n=") + std::to_string(n)));
}

std::vector<float> regSpace::axisToKey(const int ix) const {
  int iaxis = ix - 1;
  long long naft = 1, nbef = 1, naxis;
  std::vector<axis> axes = getHyper()->getAxes();
  for (auto iax = 1; iax < axes.size(); iax++) {
    if (iax < iaxis && iax != 0)
      nbef = nbef * axes[iax].n;
    else if (iax == iaxis)
      naxis = axes[iax].n;
    else
      naft = naft * axes[iax].n;
  }

  long long i = 0;
  float o = axes[iaxis].o;
  float d = axes[iaxis].d;

  std::vector<float> key(nbef * naft * naxis);
  for (long long i3 = 0; i3 < naft; i3++) {
    for (long long i2 = 0; i2 < naxis; i2++) {
      for (long long i1 = 0; i1 < nbef; i1++, i++) {
        key[i] = o + i2 * d;
      }
    }
  }
  return key;
}
