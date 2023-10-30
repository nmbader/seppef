#include "byte1DReg.h"
#include "byte2DReg.h"
#include "byte3DReg.h"
#include "byte4DReg.h"
#include "byte5DReg.h"
#include "byte6DReg.h"
#include "complex1DReg.h"
#include "complex2DReg.h"
#include "complex3DReg.h"
#include "complex4DReg.h"
#include "complex5DReg.h"
#include "complex6DReg.h"
#include "complexDouble1DReg.h"
#include "complexDouble2DReg.h"
#include "complexDouble3DReg.h"
#include "complexDouble4DReg.h"
#include "complexDouble5DReg.h"
#include "complexDouble6DReg.h"
#include "double1DReg.h"
#include "double2DReg.h"
#include "double3DReg.h"
#include "double4DReg.h"
#include "double5DReg.h"
#include "double6DReg.h"
#include "float1DReg.h"
#include "float2DReg.h"
#include "float3DReg.h"
#include "float4DReg.h"
#include "float5DReg.h"
#include "float6DReg.h"
#include "int1DReg.h"
#include "int2DReg.h"
#include "int3DReg.h"
#include "int4DReg.h"
#include "int5DReg.h"
#include "int6DReg.h"
#include "rectFilter.h"
#include "regSpace.h"
#include "sepVectorConfig.h"
#include "short1DReg.h"
#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
namespace SEP {
using namespace SEP;

PYBIND11_MODULE(pySepVector, clsVector) {
  py::class_<Vector, std::shared_ptr<Vector>>(clsVector, "Vector")
      .def(py::init<>(), "Initlialize a new Vector (don't use this");

  py::class_<regSpace, std::shared_ptr<regSpace>>(clsVector, "regSpace")
      .def(py::init<>(), "Initlialize a new regSpace (don't use this")
      .def("getHyper",
           (std::shared_ptr<hypercube>(regSpace::*)()) & regSpace::getHyper,
           "Get the hypercube")
      .def_property("_hyper", &regSpace::getHyper, &regSpace::setHyper,
                    py::return_value_policy::reference)

      .def("setHyper",
           (void (regSpace::*)(std::shared_ptr<hypercube>)) &
               regSpace::setHyper,
           "Get the hypercube")
      .def("getVoidPtr", (void *(regSpace::*)()) & regSpace::getVoidPtr,
           "Get a void ptr")
      .def("getEsize", (int (regSpace::*)()) & regSpace::getEsize,
           "Get the element size");

  py::class_<floatHyper, Vector, regSpace, std::shared_ptr<floatHyper>>(
      clsVector,
      "floatHyper") //
      .def(py::init<>(), "Initlialize a new Float Hyper (don't use this")
      .def("getSpaceOnly", (bool (floatHyper::*)()) & floatHyper::getSpaceOnly,
           "Check to see if this only a description of the vector space")

      .def("setData", (void (floatHyper::*)(float *)) & floatHyper::setData,
           "Set the data pointer")
      .def("getVals", (float *(floatHyper::*)()) & floatHyper::getVals,
           "Get the data pointer")

      .def("isDifferent",
           (bool (floatHyper::*)(std::shared_ptr<floatHyper>)) &
               floatHyper::isDifferent,
           "Check to  see if two vectors are different")

      .def_property("_vals", &floatHyper::getVals, &floatHyper::setData,
                    py::return_value_policy::reference)

      .def("createMask",
           (void (floatHyper::*)(const float, const float)) &
               floatHyper::createMask,
           "Create a mask weight")
      .def("add",
           (void (floatHyper::*)(std::shared_ptr<floatHyper>)) &
               floatHyper::add,
           "Add two vectors")
      .def("scale", (void (floatHyper::*)(const double)) & floatHyper::scale,
           "Scale a vector")
      .def("clipVector",
           (void (floatHyper::*)(const std::shared_ptr<floatHyper>,
                                 const std::shared_ptr<floatHyper>)) &
               floatHyper::clipVector,
           "vec=min(max(low,vec),high)")
      .def("calcHisto",
           (void (floatHyper::*)(std::shared_ptr<int1DReg>, float, float)) &
               floatHyper::calcHisto,
           "Calculate histogram")
      .def("scaleAdd",
           (void (floatHyper::*)(std::shared_ptr<floatHyper>, const double,
                                 const double)) &
               floatHyper::scaleAdd,
           "vec=vec*sc1+vec2*sc2")
      .def("calcCheckSum",
           (double (floatHyper::*)() const) & floatHyper::calcCheckSum,
           "Calculate checksum of a vector")
      .def("clip",
           (void (floatHyper::*)(const float, const float)) & floatHyper::clip,
           "Clip a dataset given minumum and maximum value")
      .def("cent",
           (float (floatHyper::*)(const float, const int) const) &
               floatHyper::cent,
           "Calculate  the percentile of a dataset")
      .def("norm",
           (double (floatHyper::*)(const int n) const) & floatHyper::norm,
           "Calculate n-norm of a vector")
      .def("zero", (void (floatHyper::*)()) & floatHyper::zero,
           "Fill a vector with zero")
      .def("set", (void (floatHyper::*)(const float)) & floatHyper::set,
           "Fill a vector with a value")
      .def("signum", (void (floatHyper::*)()) & floatHyper::signum,
           "Fill a vector with it's signum")
      .def("rand", (void (floatHyper::*)()) & floatHyper::random,
           "Fill a vector with random number")

      .def("mult",
           (void (floatHyper::*)(std::shared_ptr<floatHyper>)) &
               floatHyper::mult,
           "vec=vec*vec2")
      .def("dot",
           (double (floatHyper::*)(std::shared_ptr<floatHyper>) const) &
               floatHyper::dot,
           "Calculate dot product")
      .def("checkSame",
           (bool (floatHyper::*)(const std::shared_ptr<floatHyper>) const) &
               floatHyper::checkSame,
           "Check to make sure the vectors exist in the same space");

  py::class_<float1DReg, floatHyper, std::shared_ptr<float1DReg>>(
      clsVector, "float1DReg", py::buffer_protocol())
      .def(py::init<const int>(), "Initialize giving size")
      .def(py::init<const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const std::shared_ptr<float6DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<float5DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<float4DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<float3DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<float2DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 2-D hypercube")
      .def(py::init<const std::shared_ptr<float1DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 1-D hypercube")
      .def("clone",
           (std::shared_ptr<float1DReg>(float1DReg::*)() const) &
               float1DReg::clone,
           "Make a copy of the vector")
      .def("window",
           (std::shared_ptr<float1DReg>(float1DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               float1DReg::window,
           "Window a vector")

      .def("cloneSpace",
           (std::shared_ptr<float1DReg>(float1DReg::*)() const) &
               float1DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (float1DReg::*)()) & float1DReg::allocate,
           "Allocate the array")

      .def_buffer([](float1DReg &m) -> py::buffer_info {
        return py::buffer_info(m.getVals(), sizeof(float),
                               py::format_descriptor<float>::format(), 1,
                               {m.getHyper()->getAxis(1).n}, {sizeof(float)});
      });
  py::class_<float2DReg, floatHyper, std::shared_ptr<float2DReg>>(
      clsVector, "float2DReg", py::buffer_protocol())
      .def(py::init<const int, const int>(), "Initialize giving size")
      .def(py::init<const axis &, const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const std::shared_ptr<float6DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<float5DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<float4DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<float3DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<float2DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 2-D hypercube")
      .def("clone",
           (std::shared_ptr<float2DReg>(float2DReg::*)() const) &
               float2DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<float2DReg>(float2DReg::*)() const) &
               float2DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (float2DReg::*)()) & float2DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<float2DReg>(float2DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               float2DReg::window,
           "Window a vector")

      .def_buffer([](float2DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(float), py::format_descriptor<float>::format(),
            2, {m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(float) * m.getHyper()->getAxis(1).n, sizeof(float)});
      });

  py::class_<float3DReg, floatHyper, std::shared_ptr<float3DReg>>(
      clsVector, "float3DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def("allocate", (void (float3DReg::*)()) & float3DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<float3DReg>(float3DReg::*)() const) &
               float3DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<float3DReg>(float3DReg::*)() const) &
               float3DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<float1DReg>(float3DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               float3DReg::window,
           "Window a vector")

      .def_buffer([](float3DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(float), py::format_descriptor<float>::format(),
            3,
            {m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(float) * m.getHyper()->getAxis(1).n, sizeof(float)});
      });

  py::class_<float4DReg, floatHyper, std::shared_ptr<float4DReg>>(
      clsVector, "float4DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (float4DReg::*)()) & float4DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<float4DReg>(float4DReg::*)() const) &
               float4DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<float4DReg>(float4DReg::*)() const) &
               float4DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<float4DReg>(float4DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               float4DReg::window,
           "Window a vector")

      .def_buffer([](float4DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(float), py::format_descriptor<float>::format(),
            4,
            {m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(float) * m.getHyper()->getAxis(1).n, sizeof(float)});
      });
  py::class_<float5DReg, floatHyper, std::shared_ptr<float5DReg>>(
      clsVector, "float5DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (float5DReg::*)()) & float5DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<float5DReg>(float5DReg::*)() const) &
               float5DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<float5DReg>(float5DReg::*)() const) &
               float5DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<float5DReg>(float5DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               float5DReg::window,
           "Window a vector")

      .def_buffer([](float5DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(float), py::format_descriptor<float>::format(),
            5,
            {m.getHyper()->getAxis(5).n, m.getHyper()->getAxis(4).n,
             m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(float) * m.getHyper()->getAxis(1).n, sizeof(float)});
      });

  py::class_<float6DReg, floatHyper, std::shared_ptr<float6DReg>>(
      clsVector, "float6DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int,
                    const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (float6DReg::*)()) & float6DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<float6DReg>(float6DReg::*)() const) &
               float6DReg::clone,
           "Make a copy of the vector")
      .def("window",
           (std::shared_ptr<float6DReg>(float6DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               float6DReg::window,
           "Window a vector")

      .def("cloneSpace",
           (std::shared_ptr<float6DReg>(float6DReg::*)() const) &
               float6DReg::cloneSpace,
           "Make a copy of the vector space")
      .def_buffer([](float6DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(float), py::format_descriptor<float>::format(),
            6,
            {m.getHyper()->getAxis(6).n, m.getHyper()->getAxis(5).n,
             m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n * m.getHyper()->getAxis(5).n,
             sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(float) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(float) * m.getHyper()->getAxis(1).n, sizeof(float)});
      });

  py::class_<rectFilter2D, float2DReg, std::shared_ptr<rectFilter2D>>(
      clsVector, "rectFilter2D")
      .def(py::init<const std::vector<int> &, const std::vector<int> &,
                    const bool>(),
           "Initialize rectFilter2D")
      .def("clone",
           (std::shared_ptr<rectFilter2D>(rectFilter2D::*)() const) &
               rectFilter2D::clone,
           "Make a copy of the vector")

      .def("cloneSpace",
           (std::shared_ptr<rectFilter2D>(rectFilter2D::*)() const) &
               rectFilter2D::cloneSpace,
           "Make a copy of the vector space")
      .def_buffer([](rectFilter2D &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(float), py::format_descriptor<float>::format(),
            2,
            {
                m.getHyper()->getAxis(2).n,
                m.getHyper()->getAxis(1).n,
            },
            {sizeof(float) * m.getHyper()->getAxis(1).n, sizeof(float)});
      });

  py::class_<rectFilter1D, float1DReg, std::shared_ptr<rectFilter1D>>(
      clsVector, "rectFilter1D", py::buffer_protocol())
      .def(py::init<const std::vector<int> &, const std::vector<int> &,
                    const bool>(),
           "Initialize rectFilter1D")
      .def("clone",
           (std::shared_ptr<rectFilter1D>(rectFilter1D::*)() const) &
               rectFilter1D::clone,
           "Make a copy of the vector")

      .def("cloneSpace",
           (std::shared_ptr<rectFilter1D>(rectFilter1D::*)() const) &
               rectFilter1D::cloneSpace,
           "Make a copy of the vector space")
      .def_buffer([](rectFilter1D &m) -> py::buffer_info {
        return py::buffer_info(m.getVals(), sizeof(float),
                               py::format_descriptor<float>::format(), 1,
                               {m.getHyper()->getAxis(1).n}, {sizeof(float)});
      });

  py::class_<doubleHyper, regSpace, std::shared_ptr<doubleHyper>>(
      clsVector,
      "doubleHyper") //
      .def(py::init<>(), "Initlialize a new double Hyper (don't use this")
      .def("getSpaceOnly",
           (bool (doubleHyper::*)()) & doubleHyper::getSpaceOnly,
           "Check to see if this only a description of the vector space")

      .def("setData", (void (doubleHyper::*)(double *)) & doubleHyper::setData,
           "Set the data pointer")
      .def("getVals", (double *(doubleHyper::*)()) & doubleHyper::getVals,
           "Get the data pointer")
      .def("clip",
           (void (doubleHyper::*)(const double, const double)) &
               doubleHyper::clip,
           "Clip a dataset given minumum and maximum value")
      .def("cent",
           (double (doubleHyper::*)(const float pct, const int jsamp) const) &
               doubleHyper::cent,
           "Calculate  the percentile of a dataset")
      .def("calcHisto",
           (void (doubleHyper::*)(std::shared_ptr<int1DReg>, float, float)) &
               doubleHyper::calcHisto,
           "Calculate histogram")
      .def("clipVector",
           (void (doubleHyper::*)(const std::shared_ptr<doubleHyper>,
                                  const std::shared_ptr<doubleHyper>)) &
               doubleHyper::clipVector,
           "vec=min(max(low,vec),high)")
      .def("isDifferent",
           (bool (doubleHyper::*)(std::shared_ptr<doubleHyper>)) &
               doubleHyper::isDifferent,
           "Check to  see if two vectors are different")

      .def_property("_vals", &doubleHyper::getVals, &doubleHyper::setData,
                    py::return_value_policy::reference)

      .def("createMask",
           (void (doubleHyper::*)(const double, const double)) &
               doubleHyper::createMask,
           "Create a mask weight")
      .def("add",
           (void (doubleHyper::*)(std::shared_ptr<doubleHyper>)) &
               doubleHyper::add,
           "Add two vectors")
      .def("scale", (void (doubleHyper::*)(const double)) & doubleHyper::scale,
           "Scale a vector")

      .def("scaleAdd",
           (void (doubleHyper::*)(std::shared_ptr<doubleHyper>, const double,
                                  const double)) &
               doubleHyper::scaleAdd,
           "vec=vec*sc1+vec2*sc2")
      .def("calcCheckSum",
           (double (doubleHyper::*)() const) & doubleHyper::calcCheckSum,
           "Calculate checksum of a vector")

      .def("norm",
           (double (doubleHyper::*)(const int n) const) & doubleHyper::norm,
           "Calculate n-norm of a vector")
      .def("zero", (void (doubleHyper::*)()) & doubleHyper::zero,
           "Fill a vector with zero")
      .def("set", (void (doubleHyper::*)(const double)) & doubleHyper::set,
           "Fill a vector with a value")
      .def("signum", (void (doubleHyper::*)()) & doubleHyper::signum,
           "Fill a vector with it's signum")
      .def("rand", (void (doubleHyper::*)()) & doubleHyper::random,
           "Fill a vector with random number")

      .def("mult",
           (void (doubleHyper::*)(std::shared_ptr<doubleHyper>)) &
               doubleHyper::mult,
           "vec=vec*vec2")
      .def("dot",
           (double (doubleHyper::*)(std::shared_ptr<doubleHyper>) const) &
               doubleHyper::dot,
           "Calculate dot product")
      .def("checkSame",
           (bool (doubleHyper::*)(const std::shared_ptr<doubleHyper>) const) &
               doubleHyper::checkSame,
           "Check to make sure the vectors exist in the same space");

  py::class_<double1DReg, doubleHyper, std::shared_ptr<double1DReg>>(
      clsVector, "double1DReg", py::buffer_protocol())
      .def(py::init<const int>(), "Initialize giving size")
      .def(py::init<const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const std::shared_ptr<double6DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<double5DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<double4DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<double3DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<double2DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 2-D hypercube")
      .def(py::init<const std::shared_ptr<double1DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 1-D hypercube")
      .def("clone",
           (std::shared_ptr<double1DReg>(double1DReg::*)() const) &
               double1DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<double1DReg>(double1DReg::*)() const) &
               double1DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (double1DReg::*)()) & double1DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<double1DReg>(double1DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               double1DReg::window,
           "Window a vector")

      .def_buffer([](double1DReg &m) -> py::buffer_info {
        return py::buffer_info(m.getVals(), sizeof(double),
                               py::format_descriptor<double>::format(), 1,
                               {m.getHyper()->getAxis(1).n}, {sizeof(double)});
      });
  py::class_<double2DReg, doubleHyper, std::shared_ptr<double2DReg>>(
      clsVector, "double2DReg", py::buffer_protocol())
      .def(py::init<const int, const int>(), "Initialize giving size")
      .def(py::init<const axis &, const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("clone",
           (std::shared_ptr<double2DReg>(double2DReg::*)() const) &
               double2DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<double2DReg>(double2DReg::*)() const) &
               double2DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (double2DReg::*)()) & double2DReg::allocate,
           "Allocate the array")
      .def(py::init<const std::shared_ptr<double6DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<double5DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<double4DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<double3DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<double2DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 2-D hypercube")
      .def("window",
           (std::shared_ptr<double2DReg>(double2DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               double2DReg::window,
           "Window a vector")

      .def_buffer([](double2DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(double),
            py::format_descriptor<double>::format(), 2,
            {m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(double) * m.getHyper()->getAxis(1).n, sizeof(double)});
      });

  py::class_<double3DReg, doubleHyper, std::shared_ptr<double3DReg>>(
      clsVector, "double3DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def("allocate", (void (double3DReg::*)()) & double3DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<double3DReg>(double3DReg::*)() const) &
               double3DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<double3DReg>(double3DReg::*)() const) &
               double3DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<double3DReg>(double3DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               double3DReg::window,
           "Window a vector")

      .def_buffer([](double3DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(double),
            py::format_descriptor<double>::format(), 3,
            {m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(double) * m.getHyper()->getAxis(1).n, sizeof(double)});
      });

  py::class_<double4DReg, doubleHyper, std::shared_ptr<double4DReg>>(
      clsVector, "double4DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (double4DReg::*)()) & double4DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<double4DReg>(double4DReg::*)() const) &
               double4DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<double4DReg>(double4DReg::*)() const) &
               double4DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<double4DReg>(double4DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               double4DReg::window,
           "Window a vector")
      .def_buffer([](double4DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(double),
            py::format_descriptor<double>::format(), 4,
            {m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(double) * m.getHyper()->getAxis(1).n, sizeof(double)});
      });
  py::class_<double5DReg, doubleHyper, std::shared_ptr<double5DReg>>(
      clsVector, "double5DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (double5DReg::*)()) & double5DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<double5DReg>(double5DReg::*)() const) &
               double5DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<double5DReg>(double5DReg::*)() const) &
               double5DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<double5DReg>(double5DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               double5DReg::window,
           "Window a vector")
      .def_buffer([](double5DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(double),
            py::format_descriptor<double>::format(), 5,
            {m.getHyper()->getAxis(5).n, m.getHyper()->getAxis(4).n,
             m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(double) * m.getHyper()->getAxis(1).n, sizeof(double)});
      });

  py::class_<double6DReg, doubleHyper, std::shared_ptr<double6DReg>>(
      clsVector, "double6DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int,
                    const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (double6DReg::*)()) & double6DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<double6DReg>(double6DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               double6DReg::window,
           "Window a vector")
      .def("clone",
           (std::shared_ptr<double6DReg>(double6DReg::*)() const) &
               double6DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<double6DReg>(double6DReg::*)() const) &
               double6DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<double6DReg>(double6DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               double6DReg::window,
           "Window a vector")

      .def_buffer([](double6DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(double),
            py::format_descriptor<double>::format(), 6,
            {m.getHyper()->getAxis(6).n, m.getHyper()->getAxis(5).n,
             m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n * m.getHyper()->getAxis(5).n,
             sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(double) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(double) * m.getHyper()->getAxis(1).n, sizeof(double)});
      });

  py::class_<shortHyper, regSpace, std::shared_ptr<shortHyper>>(clsVector,
                                                                "shortHyper") //
      .def(py::init<>(), "Initlialize a new int Hyper (don't use this")
      .def("getSpaceOnly", (bool (shortHyper::*)()) & shortHyper::getSpaceOnly,
           "Check to see if this only a description of the vector space")

      .def("setData", (void (shortHyper::*)(int *)) & shortHyper::setData,
           "Set the data pointer")
      .def("getVals", (int *(shortHyper::*)()) & shortHyper::getVals,
           "Get the data pointer")

      .def("isDifferent",
           (bool (shortHyper::*)(std::shared_ptr<shortHyper>)) &
               shortHyper::isDifferent,
           "Check to  see if two vectors are different")

      .def_property("_vals", &shortHyper::getVals, &shortHyper::setData,
                    py::return_value_policy::reference)

      .def("add",
           (void (shortHyper::*)(std::shared_ptr<shortHyper>)) &
               shortHyper::add,
           "Add two vectors")
      .def("scale", (void (shortHyper::*)(const int)) & shortHyper::scale,
           "Scale a vector")
      .def("scaleAdd",
           (void (shortHyper::*)(std::shared_ptr<shortHyper>, const int,
                                 const int)) &
               shortHyper::scaleAdd,
           "vec=vec*sc1+vec2*sc2")
      .def("calcCheckSum",
           (int (shortHyper::*)() const) & shortHyper::calcCheckSum,
           "Calculate checksum of a vector")
      .def("zero", (void (shortHyper::*)()) & shortHyper::zero,
           "Fill a vector with zero")
      .def("set", (void (shortHyper::*)(const short)) & shortHyper::set,
           "Fill a vector with a value")
      .def("signum", (void (shortHyper::*)()) & shortHyper::signum,
           "Fill a vector with it's signum")
      .def("rand", (void (shortHyper::*)()) & shortHyper::random,
           "Fill a vector with random number")
      .def("mult",
           (void (shortHyper::*)(std::shared_ptr<shortHyper>)) &
               shortHyper::mult,
           "vec=vec*vec2")
      .def("checkSame",
           (bool (shortHyper::*)(const std::shared_ptr<shortHyper>) const) &
               shortHyper::checkSame,
           "Check to make sure the vectors exist in the same space");
  py::class_<short1DReg, shortHyper, std::shared_ptr<short1DReg>>(
      clsVector, "short1DReg", py::buffer_protocol())
      .def(py::init<const int>(), "Initialize giving size")
      .def(py::init<const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")

      .def("clone",
           (std::shared_ptr<short1DReg>(short1DReg::*)() const) &
               short1DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<short1DReg>(short1DReg::*)() const) &
               short1DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (short1DReg::*)()) & short1DReg::allocate,
           "Allocate the array")

      .def_buffer([](short1DReg &m) -> py::buffer_info {
        return py::buffer_info(m.getVals(), sizeof(short),
                               py::format_descriptor<int>::format(), 1,
                               {m.getHyper()->getAxis(1).n}, {sizeof(short)});
      });

  py::class_<intHyper, regSpace, std::shared_ptr<intHyper>>(clsVector,
                                                            "intHyper") //
      .def(py::init<>(), "Initlialize a new int Hyper (don't use this")
      .def("getSpaceOnly", (bool (intHyper::*)()) & intHyper::getSpaceOnly,
           "Check to see if this only a description of the vector space")

      .def("setData", (void (intHyper::*)(int *)) & intHyper::setData,
           "Set the data pointer")
      .def("getVals", (int *(intHyper::*)()) & intHyper::getVals,
           "Get the data pointer")

      .def("isDifferent",
           (bool (intHyper::*)(std::shared_ptr<intHyper>)) &
               intHyper::isDifferent,
           "Check to  see if two vectors are different")

      .def_property("_vals", &intHyper::getVals, &intHyper::setData,
                    py::return_value_policy::reference)
      .def("clip", (void (intHyper::*)(const int, const int)) & intHyper::clip,
           "Clip a dataset given minumum and maximum value")
      .def("cent",
           (int (intHyper::*)(const float pct, const int jsamp) const) &
               intHyper::cent,
           "Calculate  the percentile of a dataset")
      .def("add",
           (void (intHyper::*)(std::shared_ptr<intHyper>)) & intHyper::add,
           "Add two vectors")
      .def("scale", (void (intHyper::*)(const int)) & intHyper::scale,
           "Scale a vector")
      .def("scaleAdd",
           (void (intHyper::*)(std::shared_ptr<intHyper>, const int,
                               const int)) &
               intHyper::scaleAdd,
           "vec=vec*sc1+vec2*sc2")
      .def("calcCheckSum", (int (intHyper::*)() const) & intHyper::calcCheckSum,
           "Calculate checksum of a vector")
      .def("zero", (void (intHyper::*)()) & intHyper::zero,
           "Fill a vector with zero")
      .def("set", (void (intHyper::*)(const int)) & intHyper::set,
           "Fill a vector with a value")
      .def("signum", (void (intHyper::*)()) & intHyper::signum,
           "Fill a vector with it's signum")
      .def("rand", (void (intHyper::*)()) & intHyper::random,
           "Fill a vector with random number")
      .def("mult",
           (void (intHyper::*)(std::shared_ptr<intHyper>)) & intHyper::mult,
           "vec=vec*vec2")
      .def("checkSame",
           (bool (intHyper::*)(const std::shared_ptr<intHyper>) const) &
               intHyper::checkSame,
           "Check to make sure the vectors exist in the same space");

  py::class_<int1DReg, intHyper, std::shared_ptr<int1DReg>>(
      clsVector, "int1DReg", py::buffer_protocol())
      .def(py::init<const int>(), "Initialize giving size")
      .def(py::init<const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const std::shared_ptr<int6DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<int5DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<int4DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<int3DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<int2DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 2-D hypercube")
      .def(py::init<const std::shared_ptr<int1DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 1-D hypercube")
      .def("clone",
           (std::shared_ptr<int1DReg>(int1DReg::*)() const) & int1DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<int1DReg>(int1DReg::*)() const) &
               int1DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (int1DReg::*)()) & int1DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<int1DReg>(int1DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               int1DReg::window,
           "Window a vector")
      .def_buffer([](int1DReg &m) -> py::buffer_info {
        return py::buffer_info(m.getVals(), sizeof(int),
                               py::format_descriptor<int>::format(), 1,
                               {m.getHyper()->getAxis(1).n}, {sizeof(int)});
      });
  py::class_<int2DReg, intHyper, std::shared_ptr<int2DReg>>(
      clsVector, "int2DReg", py::buffer_protocol())
      .def(py::init<const int, const int>(), "Initialize giving size")
      .def(py::init<const axis &, const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const std::shared_ptr<int6DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<int5DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<int4DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<int3DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<int2DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 2-D hypercube")
      .def("clone",
           (std::shared_ptr<int2DReg>(int2DReg::*)() const) & int2DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<int2DReg>(int2DReg::*)() const) &
               int2DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (int2DReg::*)()) & int2DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<int2DReg>(int2DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               int2DReg::window,
           "Window a vector")

      .def_buffer([](int2DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(int), py::format_descriptor<int>::format(), 2,
            {m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(int) * m.getHyper()->getAxis(1).n, sizeof(int)});
      });

  py::class_<int3DReg, intHyper, std::shared_ptr<int3DReg>>(
      clsVector, "int3DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def("allocate", (void (int3DReg::*)()) & int3DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<int3DReg>(int3DReg::*)() const) & int3DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<int3DReg>(int3DReg::*)() const) &
               int3DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<int3DReg>(int3DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               int3DReg::window,
           "Window a vector")

      .def_buffer([](int3DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(int), py::format_descriptor<int>::format(), 3,
            {m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(int) * m.getHyper()->getAxis(1).n, sizeof(int)});
      });

  py::class_<int4DReg, intHyper, std::shared_ptr<int4DReg>>(
      clsVector, "int4DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (int4DReg::*)()) & int4DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<int4DReg>(int4DReg::*)() const) & int4DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<int4DReg>(int4DReg::*)() const) &
               int4DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<int4DReg>(int4DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               int4DReg::window,
           "Window a vector")
      .def_buffer([](int4DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(int), py::format_descriptor<int>::format(), 4,
            {m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(int) * m.getHyper()->getAxis(1).n, sizeof(int)});
      });
  py::class_<int5DReg, intHyper, std::shared_ptr<int5DReg>>(
      clsVector, "int5DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (int5DReg::*)()) & int5DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<int5DReg>(int5DReg::*)() const) & int5DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<int5DReg>(int5DReg::*)() const) &
               int5DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<int5DReg>(int5DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               int5DReg::window,
           "Window a vector")
      .def_buffer([](int5DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(int), py::format_descriptor<int>::format(), 5,
            {m.getHyper()->getAxis(5).n, m.getHyper()->getAxis(4).n,
             m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(int) * m.getHyper()->getAxis(1).n, sizeof(int)});
      });

  py::class_<int6DReg, intHyper, std::shared_ptr<int6DReg>>(
      clsVector, "int6DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int,
                    const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (int6DReg::*)()) & int6DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<int6DReg>(int6DReg::*)() const) & int6DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<int6DReg>(int6DReg::*)() const) &
               int6DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<int6DReg>(int6DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               int6DReg::window,
           "Window a vector")
      .def_buffer([](int6DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(int), py::format_descriptor<int>::format(), 6,
            {m.getHyper()->getAxis(6).n, m.getHyper()->getAxis(5).n,
             m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n * m.getHyper()->getAxis(5).n,
             sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(int) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(int) * m.getHyper()->getAxis(1).n, sizeof(int)});
      });

  py::class_<byteHyper, regSpace, std::shared_ptr<byteHyper>>(clsVector,
                                                              "byteHyper") //
      .def(py::init<>(), "Initlialize a new byte Hyper (don't use this")
      .def("getSpaceOnly", (bool (byteHyper::*)()) & byteHyper::getSpaceOnly,
           "Check to see if this only a description of the vector space")

      .def("setData",
           (void (byteHyper::*)(unsigned char *)) & byteHyper::setData,
           "Set the data pointer")
      .def("getVals", (unsigned char *(byteHyper::*)()) & byteHyper::getVals,
           "Get the data pointer")
      .def("calcHisto",
           (void (byteHyper::*)(std::shared_ptr<int1DReg>, float, float)) &
               byteHyper::calcHisto,
           "Calculate histogram")
      .def("isDifferent",
           (bool (byteHyper::*)(std::shared_ptr<byteHyper>)) &
               byteHyper::isDifferent,
           "Check to  see if two vectors are different")

      .def_property("_vals", &byteHyper::getVals, &byteHyper::setData,
                    py::return_value_policy::reference)
      .def("clip",
           (void (byteHyper::*)(const int, const int)) & byteHyper::clip,
           "Clip a dataset given minumum and maximum value")
      .def("cent",
           (unsigned char (byteHyper::*)(const float, const int) const) &
               byteHyper::cent,
           "Calculate  the percentile of a dataset")
      .def("calcCheckSum",
           (unsigned char (byteHyper::*)() const) & byteHyper::calcCheckSum,
           "Calculate checksum of a vector")

      .def("zero", (void (byteHyper::*)()) & byteHyper::zero,
           "Fill a vector with zero")
      .def("set", (void (byteHyper::*)(const unsigned char)) & byteHyper::set,
           "Fill a vector with a value")
      .def("rand", (void (byteHyper::*)()) & byteHyper::random,
           "Fill a vector with random number")

      .def("checkSame",
           (bool (byteHyper::*)(const std::shared_ptr<byteHyper>) const) &
               byteHyper::checkSame,
           "Check to make sure the vectors exist in the same space");

  py::class_<byte1DReg, byteHyper, std::shared_ptr<byte1DReg>>(
      clsVector, "byte1DReg", py::buffer_protocol())
      .def(py::init<const int>(), "Initialize giving size")
      .def(py::init<const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const std::shared_ptr<byte6DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<byte5DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<byte4DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<byte3DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<byte2DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 2-D hypercube")
      .def(py::init<const std::shared_ptr<byte1DReg>, const int, const bool,
                    const std::vector<int> &, const std::vector<int> &,
                    std::vector<int> &>(),
           "Create a 1-D slice from 1-D hypercube")
      .def("clone",
           (std::shared_ptr<byte1DReg>(byte1DReg::*)() const) &
               byte1DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<byte1DReg>(byte1DReg::*)() const) &
               byte1DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (byte1DReg::*)()) & byte1DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<byte1DReg>(byte1DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               byte1DReg::window,
           "Window a vector")

      .def_buffer([](byte1DReg &m) -> py::buffer_info {
        return py::buffer_info(m.getVals(), sizeof(unsigned char),
                               py::format_descriptor<unsigned char>::format(),
                               1, {m.getHyper()->getAxis(1).n},
                               {sizeof(unsigned char)});
      });
  py::class_<byte2DReg, byteHyper, std::shared_ptr<byte2DReg>>(
      clsVector, "byte2DReg", py::buffer_protocol())
      .def(py::init<const int, const int>(), "Initialize giving size")
      .def(py::init<const axis &, const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const std::shared_ptr<byte6DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 6-D hypercube")
      .def(py::init<const std::shared_ptr<byte5DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 5-D hypercube")
      .def(py::init<const std::shared_ptr<byte4DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 4-D hypercube")
      .def(py::init<const std::shared_ptr<byte3DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 3-D hypercube")
      .def(py::init<const std::shared_ptr<byte2DReg>, const int, const bool,
                    const int, const bool, const std::vector<int> &,
                    const std::vector<int> &, std::vector<int> &>(),
           "Create a 2-D slice from 2-D hypercube")
      .def("clone",
           (std::shared_ptr<byte2DReg>(byte2DReg::*)() const) &
               byte2DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<byte2DReg>(byte2DReg::*)() const) &
               byte2DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (byte2DReg::*)()) & byte2DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<byte2DReg>(byte2DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               byte2DReg::window,
           "Window a vector")
      .def_buffer([](byte2DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(unsigned char),
            py::format_descriptor<unsigned char>::format(), 2,
            {m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(unsigned char) * m.getHyper()->getAxis(1).n,
             sizeof(unsigned char)});
      });

  py::class_<byte3DReg, byteHyper, std::shared_ptr<byte3DReg>>(
      clsVector, "byte3DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def("allocate", (void (byte3DReg::*)()) & byte3DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<byte3DReg>(byte3DReg::*)() const) &
               byte3DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<byte3DReg>(byte3DReg::*)() const) &
               byte3DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<byte3DReg>(byte3DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               byte3DReg::window,
           "Window a vector")
      .def_buffer([](byte3DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(unsigned char),
            py::format_descriptor<unsigned char>::format(), 3,
            {m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n,
             sizeof(unsigned char)});
      });

  py::class_<byte4DReg, byteHyper, std::shared_ptr<byte4DReg>>(
      clsVector, "byte4DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (byte4DReg::*)()) & byte4DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<byte4DReg>(byte4DReg::*)() const) &
               byte4DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<byte4DReg>(byte4DReg::*)() const) &
               byte4DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<byte4DReg>(byte4DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               byte4DReg::window,
           "Window a vector")
      .def_buffer([](byte4DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(unsigned char),
            py::format_descriptor<unsigned char>::format(), 4,
            {m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n,
             sizeof(unsigned char)});
      });
  py::class_<byte5DReg, byteHyper, std::shared_ptr<byte5DReg>>(
      clsVector, "byte5DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (byte5DReg::*)()) & byte5DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<byte5DReg>(byte5DReg::*)() const) &
               byte5DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<byte5DReg>(byte5DReg::*)() const) &
               byte5DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<byte5DReg>(byte5DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               byte5DReg::window,
           "Window a vector")
      .def_buffer([](byte5DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(unsigned char),
            py::format_descriptor<unsigned char>::format(), 5,
            {m.getHyper()->getAxis(5).n, m.getHyper()->getAxis(4).n,
             m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n,
             sizeof(unsigned char)});
      });

  py::class_<byte6DReg, byteHyper, std::shared_ptr<byte6DReg>>(
      clsVector, "byte6DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int,
                    const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (byte6DReg::*)()) & byte6DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<byte6DReg>(byte6DReg::*)() const) &
               byte6DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<byte6DReg>(byte6DReg::*)() const) &
               byte6DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<byte6DReg>(byte6DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               byte6DReg::window,
           "Window a vector")
      .def_buffer([](byte6DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(unsigned char),
            py::format_descriptor<unsigned char>::format(), 6,
            {m.getHyper()->getAxis(6).n, m.getHyper()->getAxis(5).n,
             m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n * m.getHyper()->getAxis(5).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(unsigned char) * m.getHyper()->getAxis(1).n,
             sizeof(unsigned char)});
      });

  py::class_<complexHyper, regSpace, std::shared_ptr<complexHyper>>(
      clsVector,
      "complexHyper") //
      .def(py::init<>(), "Initlialize a new complex Hyper (don't use this")
      .def("getSpaceOnly",
           (bool (complexHyper::*)()) & complexHyper::getSpaceOnly,
           "Check to see if this only a description of the vector space")

      .def("setData",
           (void (complexHyper::*)(std::complex<float> *)) &
               complexHyper::setData,
           "Set the data pointer")
      .def("getVals",
           (std::complex<float> * (complexHyper::*)()) & complexHyper::getVals,
           "Get the data pointer")

      .def("isDifferent",
           (bool (complexHyper::*)(std::shared_ptr<complexHyper>)) &
               complexHyper::isDifferent,
           "Check to  see if two vectors are different")

      .def("mult",
           (void (complexHyper::*)(std::shared_ptr<complexHyper>)) &
               complexHyper::mult,
           "vec=vec*vec2")
      .def("dot",
           (double (complexHyper::*)(std::shared_ptr<complexHyper>) const) &
               complexHyper::dot,
           "Calculate dot product")
      .def_property("_vals", &complexHyper::getVals, &complexHyper::setData,
                    py::return_value_policy::reference)
      .def("scale",
           (void (complexHyper::*)(const double)) & complexHyper::scale,
           "Scale a vector")
      .def("scaleAdd",
           (void (complexHyper::*)(std::shared_ptr<complexHyper>, const double,
                                   const double)) &
               complexHyper::scaleAdd,
           "vec=vec*sc1+vec2*sc2")
      .def("calcCheckSum",
           (unsigned char (complexHyper::*)() const) &
               complexHyper::calcCheckSum,
           "Calculate checksum of a vector")
      .def("norm",
           (double (complexHyper::*)(const int n) const) & complexHyper::norm,
           "Calculate n-norm of a vector")
      .def("zero", (void (complexHyper::*)()) & complexHyper::zero,
           "Fill a vector with zero")
      .def("set",
           (void (complexHyper::*)(const std::complex<float>)) &
               complexHyper::set,
           "Fill a vector with a value")
      .def("rand", (void (complexHyper::*)()) & complexHyper::random,
           "Fill a vector with random number")
      .def("clipVector",
           (void (complexHyper::*)(const std::shared_ptr<floatHyper>,
                                   const std::shared_ptr<floatHyper>)) &
               complexHyper::clipVector,
           "vec=min(max(low,vec),high)")
      .def("clipVector",
           (void (complexHyper::*)(const std::shared_ptr<complexHyper>,
                                   const std::shared_ptr<complexHyper>)) &
               complexHyper::clipVector,
           "vec=min(max(low,vec),high)")
      .def("checkSame",
           (bool (complexHyper::*)(const std::shared_ptr<complexHyper>) const) &
               complexHyper::checkSame,
           "Check to make sure the vectors exist in the same space");

  py::class_<complex1DReg, complexHyper, std::shared_ptr<complex1DReg>>(
      clsVector, "complex1DReg", py::buffer_protocol())
      .def(py::init<const int>(), "Initialize giving size")
      .def(py::init<const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")

      .def("clone",
           (std::shared_ptr<complex1DReg>(complex1DReg::*)() const) &
               complex1DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<complex1DReg>(complex1DReg::*)() const) &
               complex1DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (complex1DReg::*)()) & complex1DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<complex1DReg>(complex1DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complex1DReg::window,
           "Window a vector")

      .def_buffer([](complex1DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<float>),
            py::format_descriptor<std::complex<float>>::format(), 1,
            {m.getHyper()->getAxis(1).n}, {sizeof(std::complex<float>)});
      });
  py::class_<complex2DReg, complexHyper, std::shared_ptr<complex2DReg>>(
      clsVector, "complex2DReg", py::buffer_protocol())
      .def(py::init<const int, const int>(), "Initialize giving size")
      .def(py::init<const axis &, const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("clone",
           (std::shared_ptr<complex2DReg>(complex2DReg::*)() const) &
               complex2DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<complex2DReg>(complex2DReg::*)() const) &
               complex2DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("allocate", (void (complex2DReg::*)()) & complex2DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<complex2DReg>(complex2DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complex2DReg::window,
           "Window a vector")
      .def_buffer([](complex2DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<float>),
            py::format_descriptor<std::complex<float>>::format(), 2,
            {m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<float>)});
      });

  py::class_<complex3DReg, complexHyper, std::shared_ptr<complex3DReg>>(
      clsVector, "complex3DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def("allocate", (void (complex3DReg::*)()) & complex3DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<complex3DReg>(complex3DReg::*)() const) &
               complex3DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<complex3DReg>(complex3DReg::*)() const) &
               complex3DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complex3DReg>(complex3DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complex3DReg::window,
           "Window a vector")
      .def_buffer([](complex3DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<float>),
            py::format_descriptor<std::complex<float>>::format(), 3,
            {m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<float>)});
      });

  py::class_<complex4DReg, complexHyper, std::shared_ptr<complex4DReg>>(
      clsVector, "complex4DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (complex4DReg::*)()) & complex4DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<complex4DReg>(complex4DReg::*)() const) &
               complex4DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<complex4DReg>(complex4DReg::*)() const) &
               complex4DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complex4DReg>(complex4DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complex4DReg::window,
           "Window a vector")
      .def_buffer([](complex4DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<float>),
            py::format_descriptor<std::complex<float>>::format(), 4,
            {m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<float>)});
      });
  py::class_<complex5DReg, complexHyper, std::shared_ptr<complex5DReg>>(
      clsVector, "complex5DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (complex5DReg::*)()) & complex5DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<complex5DReg>(complex5DReg::*)() const) &
               complex5DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<complex5DReg>(complex5DReg::*)() const) &
               complex5DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complex5DReg>(complex5DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complex5DReg::window,
           "Window a vector")
      .def_buffer([](complex5DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<float>),
            py::format_descriptor<std::complex<float>>::format(), 5,
            {m.getHyper()->getAxis(5).n, m.getHyper()->getAxis(4).n,
             m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<float>)});
      });

  py::class_<complex6DReg, complexHyper, std::shared_ptr<complex6DReg>>(
      clsVector, "complex6DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int,
                    const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate", (void (complex6DReg::*)()) & complex6DReg::allocate,
           "Allocate the array")
      .def("clone",
           (std::shared_ptr<complex6DReg>(complex6DReg::*)() const) &
               complex6DReg::clone,
           "Make a copy of the vector")
      .def("cloneSpace",
           (std::shared_ptr<complex6DReg>(complex6DReg::*)() const) &
               complex6DReg::cloneSpace,
           "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complex6DReg>(complex6DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complex6DReg::window,
           "Window a vector")
      .def_buffer([](complex6DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<float>),
            py::format_descriptor<std::complex<float>>::format(), 6,
            {m.getHyper()->getAxis(6).n, m.getHyper()->getAxis(5).n,
             m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n * m.getHyper()->getAxis(5).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<float>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<float>)});
      });
  py::class_<complexDoubleHyper, regSpace, std::shared_ptr<complexDoubleHyper>>(
      clsVector,
      "complexDoubleHyper") //
      .def(py::init<>(), "Initlialize a new complex Hyper (don't use this")
      .def("getSpaceOnly",
           (bool (complexDoubleHyper::*)()) & complexDoubleHyper::getSpaceOnly,
           "Check to see if this only a description of the vector space")

      .def("setData",
           (void (complexDoubleHyper::*)(std::complex<double> *)) &
               complexDoubleHyper::setData,
           "Set the data pointer")
      .def("getVals",
           (std::complex<double> * (complexDoubleHyper::*)()) &
               complexDoubleHyper::getVals,
           "Get the data pointer")

      .def("isDifferent",
           (bool (complexDoubleHyper::*)(std::shared_ptr<complexDoubleHyper>)) &
               complexDoubleHyper::isDifferent,
           "Check to  see if two vectors are different")

      .def("mult",
           (void (complexDoubleHyper::*)(std::shared_ptr<complexDoubleHyper>)) &
               complexDoubleHyper::mult,
           "vec=vec*vec2")
      .def("dot",
           (double (complexDoubleHyper::*)(std::shared_ptr<complexDoubleHyper>)
                const) &
               complexDoubleHyper::dot,
           "Calculate dot product")
      .def_property("_vals", &complexDoubleHyper::getVals,
                    &complexDoubleHyper::setData,
                    py::return_value_policy::reference)
      .def("scale",
           (void (complexDoubleHyper::*)(const double)) &
               complexDoubleHyper::scale,
           "Scale a vector")
      .def("scaleAdd",
           (void (complexDoubleHyper::*)(std::shared_ptr<complexDoubleHyper>,
                                         const double, const double)) &
               complexDoubleHyper::scaleAdd,
           "vec=vec*sc1+vec2*sc2")
      .def("calcCheckSum",
           (unsigned char (complexDoubleHyper::*)() const) &
               complexDoubleHyper::calcCheckSum,
           "Calculate checksum of a vector")
      .def("norm",
           (double (complexDoubleHyper::*)(const int n) const) &
               complexDoubleHyper::norm,
           "Calculate n-norm of a vector")
      .def("zero", (void (complexDoubleHyper::*)()) & complexDoubleHyper::zero,
           "Fill a vector with zero")
      .def("set",
           (void (complexDoubleHyper::*)(const std::complex<double>)) &
               complexDoubleHyper::set,
           "Fill a vector with a value")
      .def("rand",
           (void (complexDoubleHyper::*)()) & complexDoubleHyper::random,
           "Fill a vector with random number")
      .def("clipVector",
           (void (complexDoubleHyper::*)(const std::shared_ptr<doubleHyper>,
                                         const std::shared_ptr<doubleHyper>)) &
               complexDoubleHyper::clipVector,
           "vec=min(max(low,vec),high)")
      .def("clipVector",
           (void (complexDoubleHyper::*)(
               const std::shared_ptr<complexDoubleHyper>,
               const std::shared_ptr<complexDoubleHyper>)) &
               complexDoubleHyper::clipVector,
           "vec=min(max(low,vec),high)")
      .def("checkSame",
           (bool (complexDoubleHyper::*)(
               const std::shared_ptr<complexDoubleHyper>) const) &
               complexDoubleHyper::checkSame,
           "Check to make sure the vectors exist in the same space");

  py::class_<complexDouble1DReg, complexDoubleHyper,
             std::shared_ptr<complexDouble1DReg>>(
      clsVector, "complexDouble1DReg", py::buffer_protocol())
      .def(py::init<const int>(), "Initialize giving size")
      .def(py::init<const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")

      .def(
          "clone",
          (std::shared_ptr<complexDouble1DReg>(complexDouble1DReg::*)() const) &
              complexDouble1DReg::clone,
          "Make a copy of the vector")
      .def(
          "cloneSpace",
          (std::shared_ptr<complexDouble1DReg>(complexDouble1DReg::*)() const) &
              complexDouble1DReg::cloneSpace,
          "Make a copy of the vector space")
      .def("allocate",
           (void (complexDouble1DReg::*)()) & complexDouble1DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<complexDouble1DReg>(complexDouble1DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complexDouble1DReg::window,
           "Window a vector")

      .def_buffer([](complexDouble1DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<double>),
            py::format_descriptor<std::complex<double>>::format(), 1,
            {m.getHyper()->getAxis(1).n}, {sizeof(std::complex<double>)});
      });
  py::class_<complexDouble2DReg, complexDoubleHyper,
             std::shared_ptr<complexDouble2DReg>>(
      clsVector, "complexDouble2DReg", py::buffer_protocol())
      .def(py::init<const int, const int>(), "Initialize giving size")
      .def(py::init<const axis &, const axis &>(), "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(
          "clone",
          (std::shared_ptr<complexDouble2DReg>(complexDouble2DReg::*)() const) &
              complexDouble2DReg::clone,
          "Make a copy of the vector")
      .def(
          "cloneSpace",
          (std::shared_ptr<complexDouble2DReg>(complexDouble2DReg::*)() const) &
              complexDouble2DReg::cloneSpace,
          "Make a copy of the vector space")
      .def("allocate",
           (void (complexDouble2DReg::*)()) & complexDouble2DReg::allocate,
           "Allocate the array")
      .def("window",
           (std::shared_ptr<complexDouble2DReg>(complexDouble2DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complexDouble2DReg::window,
           "Window a vector")
      .def_buffer([](complexDouble2DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<double>),
            py::format_descriptor<std::complex<double>>::format(), 2,
            {m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<double>)});
      });

  py::class_<complexDouble3DReg, complexDoubleHyper,
             std::shared_ptr<complexDouble3DReg>>(
      clsVector, "complexDouble3DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def(py::init<const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def("allocate",
           (void (complexDouble3DReg::*)()) & complexDouble3DReg::allocate,
           "Allocate the array")
      .def(
          "clone",
          (std::shared_ptr<complexDouble3DReg>(complexDouble3DReg::*)() const) &
              complexDouble3DReg::clone,
          "Make a copy of the vector")
      .def(
          "cloneSpace",
          (std::shared_ptr<complexDouble3DReg>(complexDouble3DReg::*)() const) &
              complexDouble3DReg::cloneSpace,
          "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complexDouble3DReg>(complexDouble3DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complexDouble3DReg::window,
           "Window a vector")
      .def_buffer([](complexDouble3DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<double>),
            py::format_descriptor<std::complex<double>>::format(), 3,
            {m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<double>)});
      });

  py::class_<complexDouble4DReg, complexDoubleHyper,
             std::shared_ptr<complexDouble4DReg>>(
      clsVector, "complexDouble4DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate",
           (void (complexDouble4DReg::*)()) & complexDouble4DReg::allocate,
           "Allocate the array")
      .def(
          "clone",
          (std::shared_ptr<complexDouble4DReg>(complexDouble4DReg::*)() const) &
              complexDouble4DReg::clone,
          "Make a copy of the vector")
      .def(
          "cloneSpace",
          (std::shared_ptr<complexDouble4DReg>(complexDouble4DReg::*)() const) &
              complexDouble4DReg::cloneSpace,
          "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complexDouble4DReg>(complexDouble4DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complexDouble4DReg::window,
           "Window a vector")
      .def_buffer([](complexDouble4DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<double>),
            py::format_descriptor<std::complex<double>>::format(), 4,
            {m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<double>)});
      });
  py::class_<complexDouble5DReg, complexDoubleHyper,
             std::shared_ptr<complexDouble5DReg>>(
      clsVector, "complexDouble5DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate",
           (void (complexDouble5DReg::*)()) & complexDouble5DReg::allocate,
           "Allocate the array")
      .def(
          "clone",
          (std::shared_ptr<complexDouble5DReg>(complexDouble5DReg::*)() const) &
              complexDouble5DReg::clone,
          "Make a copy of the vector")
      .def(
          "cloneSpace",
          (std::shared_ptr<complexDouble5DReg>(complexDouble5DReg::*)() const) &
              complexDouble5DReg::cloneSpace,
          "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complexDouble5DReg>(complexDouble5DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complexDouble5DReg::window,
           "Window a vector")
      .def_buffer([](complexDouble5DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<double>),
            py::format_descriptor<std::complex<double>>::format(), 5,
            {m.getHyper()->getAxis(5).n, m.getHyper()->getAxis(4).n,
             m.getHyper()->getAxis(3).n, m.getHyper()->getAxis(2).n,
             m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<double>)});
      });

  py::class_<complexDouble6DReg, complexDoubleHyper,
             std::shared_ptr<complexDouble6DReg>>(
      clsVector, "complexDouble6DReg", py::buffer_protocol())
      .def(py::init<const int, const int, const int, const int, const int,
                    const int>(),
           "Initialize giving size")
      .def(py::init<const axis &, const axis &, const axis &, const axis &,
                    const axis &, const axis &>(),
           "Initialize from an axis")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize with hypercube")
      .def("allocate",
           (void (complexDouble6DReg::*)()) & complexDouble6DReg::allocate,
           "Allocate the array")
      .def(
          "clone",
          (std::shared_ptr<complexDouble6DReg>(complexDouble6DReg::*)() const) &
              complexDouble6DReg::clone,
          "Make a copy of the vector")
      .def(
          "cloneSpace",
          (std::shared_ptr<complexDouble6DReg>(complexDouble6DReg::*)() const) &
              complexDouble6DReg::cloneSpace,
          "Make a copy of the vector space")
      .def("window",
           (std::shared_ptr<complexDouble6DReg>(complexDouble6DReg::*)(
               const std::vector<int> &, const std::vector<int> &,
               const std::vector<int> &) const) &
               complexDouble6DReg::window,
           "Window a vector")
      .def_buffer([](complexDouble6DReg &m) -> py::buffer_info {
        return py::buffer_info(
            m.getVals(), sizeof(std::complex<double>),
            py::format_descriptor<std::complex<double>>::format(), 6,
            {m.getHyper()->getAxis(6).n, m.getHyper()->getAxis(5).n,
             m.getHyper()->getAxis(4).n, m.getHyper()->getAxis(3).n,
             m.getHyper()->getAxis(2).n, m.getHyper()->getAxis(1).n},
            {sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n * m.getHyper()->getAxis(5).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n *
                 m.getHyper()->getAxis(4).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n * m.getHyper()->getAxis(3).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n *
                 m.getHyper()->getAxis(2).n,
             sizeof(std::complex<double>) * m.getHyper()->getAxis(1).n,
             sizeof(std::complex<double>)});
      });
} // namespace SEP
} // namespace SEP
