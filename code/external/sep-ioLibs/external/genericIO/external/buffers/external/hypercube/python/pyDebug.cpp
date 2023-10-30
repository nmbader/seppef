#include <pybind11/pybind11.h>
#include "SEPException.h"
#include "debug.h"
namespace py = pybind11;
namespace SEP {
PYBIND11_MODULE(SEPException, m) {
  static py::exception<SEPException> exc(m, "SEPError");
  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p) std::rethrow_exception(p);
    } catch (const SEPException &e) {
      exc(e.what());
    }
  });
}
PYBIND11_MODULE(pyDebugging, m) {
  m.def("getDebug", &getDebug, "Return debugging information");
}

}  // namespace SEP
