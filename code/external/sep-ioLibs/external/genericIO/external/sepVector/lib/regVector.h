#ifndef reg_vector_h
#define reg_vector_h 1
#include <ioTypes.h>
#include "short1DReg.h"

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
namespace SEP {
/*!
   Return a vector of specifc type based on dataType and hypercube describing
   the sapce

   \param hyper Hypercube describing the space
   \param typ  Datatype (float,complex,int,..)
   \param g1   Maximum dimensions based on last axis greater 1 length
  */
std::shared_ptr<regSpace> vecFromHyper(const std::shared_ptr<hypercube> hyper,
                                       const dataType typ,
                                       const bool g1 = true);
/*!
   Return a vector that is potentially smaller than the hypercube
    \param hyper Hypercube describing larger vector space
    \param typ Datatype (float, complex, int, ...)
    \param ndim Number of dimensions for vector
    */
std::shared_ptr<regSpace> subCubeFromHyper(
    const std::shared_ptr<hypercube> hyper, const dataType typ,
    const int &ndim);
/*!
   Return a subset of the larger vector
    \param hyper Hypercube describing larger vector space
    \param nw,fw,jw Windowing parameters
    \param typ Datatype (float, complex, int, ...)
    */
std::shared_ptr<regSpace> windowFromHyper(
    const std::shared_ptr<hypercube> hyper, const std::vector<int> &nw,
    const std::vector<int> &fw, const std::vector<int> &jw, const dataType typ);
/*!
   Clone a regular space return child class
    \param storage Space to clone
    */
std::shared_ptr<regSpace> cloneRegSpace(std::shared_ptr<regSpace> storage);
}  // namespace SEP
#endif
