#include <pybind11/pybind11.h>
#include "SEPException.h"
#include "debug.h"
namespace py = pybind11;
namespace SEP {
PYBIND11_MODULE(pyDebugging, m) {
  m.def("getDebug", &getDebug, "Return debugging information");
}

}  // namespace SEP
