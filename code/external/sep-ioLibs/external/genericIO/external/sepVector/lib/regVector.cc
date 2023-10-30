#include "regVector.h"
#include "sepVectorConfig.h"
using namespace SEP;
std::shared_ptr<regSpace> SEP::vecFromHyper(
    const std::shared_ptr<hypercube> hyper, const dataType typ, const bool g1) {
  std::vector<axis> axesIn = hyper->getAxes(), axesOut;
  int ndim = hyper->getNdim();
  if (g1) ndim = hyper->getNdimG1();
  for (int i = 0; i < ndim; i++) {
    axesOut.push_back(axesIn[i]);
  }

  std::shared_ptr<hypercube> hyper2(new hypercube(axesOut));
  switch (typ) {
    case DATA_FLOAT:
      switch (ndim) {
        case 1: {
          std::shared_ptr<float1DReg> a(new float1DReg(hyper2));
          return a;
        } break;
        case 2: {
          std::shared_ptr<float2DReg> a(new float2DReg(hyper2));
          return a;
        } break;
        case 3: {
          std::shared_ptr<float3DReg> a(new float3DReg(hyper2));
          return a;
        } break;
        case 4: {
          std::shared_ptr<float4DReg> a(new float4DReg(hyper2));
          return a;
        } break;
        case 5: {
          std::shared_ptr<float5DReg> a(new float5DReg(hyper2));
          return a;
        } break;
        case 6: {
          std::shared_ptr<float6DReg> a(new float6DReg(hyper2));
          return a;
        } break;
      }
      break;
case DATA_COMPLEXDOUBLE:
      switch (ndim) {
        case 1: {
          std::shared_ptr<complexDouble1DReg> a(new complexDouble1DReg(hyper2));
          return a;
        } break;
        case 2: {
          std::shared_ptr<complexDouble2DReg> a(new complexDouble2DReg(hyper2));
          return a;
        } break;
        case 3: {
          std::shared_ptr<complexDouble3DReg> a(new complexDouble3DReg(hyper2));
          return a;
        } break;
        case 4: {
          std::shared_ptr<complexDouble4DReg> a(new complexDouble4DReg(hyper2));
          return a;
        } break;
        case 5: {
          std::shared_ptr<complexDouble5DReg> a(new complexDouble5DReg(hyper2));
          return a;
        } break;
        case 6: {
          std::shared_ptr<complexDouble6DReg> a(new complexDouble6DReg(hyper2));
          return a;
        } break;
      }

    case DATA_COMPLEX:
      switch (ndim) {
        case 1: {
          std::shared_ptr<complex1DReg> a(new complex1DReg(hyper2));
          return a;
        } break;
        case 2: {
          std::shared_ptr<complex2DReg> a(new complex2DReg(hyper2));
          return a;
        } break;
        case 3: {
          std::shared_ptr<complex3DReg> a(new complex3DReg(hyper2));
          return a;
        } break;
        case 4: {
          std::shared_ptr<complex4DReg> a(new complex4DReg(hyper2));
          return a;
        } break;
        case 5: {
          std::shared_ptr<complex5DReg> a(new complex5DReg(hyper2));
          return a;
        } break;
        case 6: {
          std::shared_ptr<complex6DReg> a(new complex6DReg(hyper2));
          return a;
        } break;
      }

      break;
    case DATA_BYTE:
      switch (ndim) {
        case 1: {
          std::shared_ptr<byte1DReg> a(new byte1DReg(hyper2));
          return a;
        } break;
        case 2: {
          std::shared_ptr<byte2DReg> a(new byte2DReg(hyper2));
          return a;
        } break;
        case 3: {
          std::shared_ptr<byte3DReg> a(new byte3DReg(hyper2));
          return a;
        } break;
        case 4: {
          std::shared_ptr<byte4DReg> a(new byte4DReg(hyper2));
          return a;
        } break;
        case 5: {
          std::shared_ptr<byte5DReg> a(new byte5DReg(hyper2));
          return a;
        } break;
        case 6: {
          std::shared_ptr<byte6DReg> a(new byte6DReg(hyper2));
          return a;
        } break;
      }

      break;
    case DATA_SHORT:
      switch (ndim) {
        case 1: {
          std::shared_ptr<short1DReg> a(new short1DReg(hyper2));
          return a;
        } break;
      }
      break;
    case DATA_INT:
      switch (ndim) {
        case 1: {
          std::shared_ptr<int1DReg> a(new int1DReg(hyper2));
          return a;
        } break;
        case 2: {
          std::shared_ptr<int2DReg> a(new int2DReg(hyper2));
          return a;
        } break;
        case 3: {
          std::shared_ptr<int3DReg> a(new int3DReg(hyper2));
          return a;
        } break;
        case 4: {
          std::shared_ptr<int4DReg> a(new int4DReg(hyper2));
          return a;
        } break;
        case 5: {
          std::shared_ptr<int5DReg> a(new int5DReg(hyper2));
          return a;
        } break;
        case 6: {
          std::shared_ptr<int6DReg> a(new int6DReg(hyper2));
          return a;
        } break;
      }

      break;
    case DATA_DOUBLE:
      switch (ndim) {
        case 1: {
          std::shared_ptr<double1DReg> a(new double1DReg(hyper2));
          return a;
        } break;
        case 2: {
          std::shared_ptr<double2DReg> a(new double2DReg(hyper2));
          return a;
        } break;
        case 3: {
          std::shared_ptr<double3DReg> a(new double3DReg(hyper2));
          return a;
        } break;
        case 4: {
          std::shared_ptr<double4DReg> a(new double4DReg(hyper2));
          return a;
        } break;
        case 5: {
          std::shared_ptr<double5DReg> a(new double5DReg(hyper2));
          return a;
        } break;
        case 6: {
          std::shared_ptr<double6DReg> a(new double6DReg(hyper2));
          return a;
        } break;
      }

      break;
    default:
      throw(SEPException(std::string("Unhandled data type")));
  }
  return nullptr;
}
std::shared_ptr<regSpace> SEP::subCubeFromHyper(
    const std::shared_ptr<hypercube> hyper, const dataType typ,
    const int &ndim) {
  std::vector<axis> axesIn = hyper->getAxes(), axesOut;

  for (int i = 0; i < ndim; i++) {
    axesOut.push_back(axesIn[i]);
  }
  std::shared_ptr<hypercube> hyper2(new hypercube(axesOut));
  return vecFromHyper(hyper2, typ);
  return nullptr;
}
/*
std::shared_ptr<regSpace> SEP::windowFromHyper(
    const std::shared_ptr<hypercube> hyper, const std::vector<int> &nw,
    const std::vector<int> &fw, const std::vector<int> &jw,
    const dataType typ) {
  size_t nbigest = 0;
  for (int i = 0; i < nw.size(); i++) {
    if (nw[i] > 1) nbigest = i;
  }
  if (hyper->getNdim() <= nbiggest)
    throw(SEPException(std::string("hyper ndim=") +
                       std::to_string(hyper->getNdim()) +
                       std::string(" nibiggest=") + std::to_string(nbiggest)));

  if (nbiggest > fw.size())
    throw(SEPException(std::string("nbiggest[") + std::to_string(nbiggest) +
                       std::string("] > fw.size()[") +
                       std::to_string(nbiggest) + std::to_string("]")));

  std::vector<axis> axes = hyper->getAxes(), axesOut;
  for (int i = 0; i < nbigest; i++) {
    checkWindow(axes[i].n, nw[i], fw[i], jw[i]);
    axis a(nw[i], axes[i].o + axes[i].d * fw[i], jw[i] * axes[i].d,
           axes[i].label);
    axesOut.push_back(a);
  }
  std::shared_ptr<hypercube> hyper2(new hypercube(axesOut));
  return vecFromHyper(hyper2, typ);
}
*/
std::shared_ptr<regSpace> SEP::cloneRegSpace(
    std::shared_ptr<regSpace> storage) {
  std::shared_ptr<float1DReg> v1 =
      std::dynamic_pointer_cast<float1DReg>(storage);
  if (v1) return v1->clone();

  std::shared_ptr<float2DReg> v2 =
      std::dynamic_pointer_cast<float2DReg>(storage);
  if (v2) return v2->clone();

  std::shared_ptr<float3DReg> v3 =
      std::dynamic_pointer_cast<float3DReg>(storage);
  if (v3) return v3->clone();

  std::shared_ptr<float4DReg> v4 =
      std::dynamic_pointer_cast<float4DReg>(storage);
  if (v4) return v4->clone();

  std::shared_ptr<float5DReg> v5 =
      std::dynamic_pointer_cast<float5DReg>(storage);
  if (v5) return v5->clone();

  std::shared_ptr<float6DReg> v6 =
      std::dynamic_pointer_cast<float6DReg>(storage);
  if (v6) return v6->clone();
  std::shared_ptr<int1DReg> va = std::dynamic_pointer_cast<int1DReg>(storage);
  if (va) return va->clone();

  std::shared_ptr<int2DReg> vb = std::dynamic_pointer_cast<int2DReg>(storage);
  if (vb) return vb->clone();

  std::shared_ptr<int3DReg> vc = std::dynamic_pointer_cast<int3DReg>(storage);
  if (vc) return vc->clone();

  std::shared_ptr<int4DReg> vd = std::dynamic_pointer_cast<int4DReg>(storage);
  if (vd) return vd->clone();

  std::shared_ptr<int5DReg> ve = std::dynamic_pointer_cast<int5DReg>(storage);
  if (ve) return ve->clone();

  std::shared_ptr<int6DReg> vf = std::dynamic_pointer_cast<int6DReg>(storage);
  if (vf) return vf->clone();

  std::shared_ptr<byte1DReg> vg = std::dynamic_pointer_cast<byte1DReg>(storage);
  if (vg) return vg->clone();

  std::shared_ptr<byte2DReg> vh = std::dynamic_pointer_cast<byte2DReg>(storage);
  if (vh) return vh->clone();

  std::shared_ptr<byte3DReg> vi = std::dynamic_pointer_cast<byte3DReg>(storage);
  if (vi) return vi->clone();

  std::shared_ptr<byte4DReg> vj = std::dynamic_pointer_cast<byte4DReg>(storage);
  if (vj) return vj->clone();

  std::shared_ptr<byte5DReg> vk = std::dynamic_pointer_cast<byte5DReg>(storage);
  if (vk) return vk->clone();

  std::shared_ptr<byte6DReg> vl = std::dynamic_pointer_cast<byte6DReg>(storage);
  if (vl) return vl->clone();

  std::shared_ptr<complex1DReg> vm =
      std::dynamic_pointer_cast<complex1DReg>(storage);
  if (vm) return vm->clone();

  std::shared_ptr<complex2DReg> vn =
      std::dynamic_pointer_cast<complex2DReg>(storage);
  if (vn) return vn->clone();

  std::shared_ptr<complex3DReg> vo =
      std::dynamic_pointer_cast<complex3DReg>(storage);
  if (vo) return vo->clone();

  std::shared_ptr<complex4DReg> vp =
      std::dynamic_pointer_cast<complex4DReg>(storage);
  if (vp) return vp->clone();

  std::shared_ptr<complex5DReg> vq =
      std::dynamic_pointer_cast<complex5DReg>(storage);
  if (vq) return vq->clone();

  std::shared_ptr<complex6DReg> vr =
      std::dynamic_pointer_cast<complex6DReg>(storage);
  if (vr) return vr->clone();
{
  std::shared_ptr<complexDouble1DReg> vm =
      std::dynamic_pointer_cast<complexDouble1DReg>(storage);
  if (vm) return vm->clone();

  std::shared_ptr<complexDouble2DReg> vn =
      std::dynamic_pointer_cast<complexDouble2DReg>(storage);
  if (vn) return vn->clone();

  std::shared_ptr<complexDouble3DReg> vo =
      std::dynamic_pointer_cast<complexDouble3DReg>(storage);
  if (vo) return vo->clone();

  std::shared_ptr<complexDouble4DReg> vp =
      std::dynamic_pointer_cast<complexDouble4DReg>(storage);
  if (vp) return vp->clone();

  std::shared_ptr<complexDouble5DReg> vq =
      std::dynamic_pointer_cast<complexDouble5DReg>(storage);
  if (vq) return vq->clone();

  std::shared_ptr<complexDouble6DReg> vr =
      std::dynamic_pointer_cast<complexDouble6DReg>(storage);
  if (vr) return vr->clone();

    }

  std::shared_ptr<double1DReg> vs =
      std::dynamic_pointer_cast<double1DReg>(storage);
  if (vs) return vs->clone();

  std::shared_ptr<double2DReg> vt =
      std::dynamic_pointer_cast<double2DReg>(storage);
  if (vt) return vt->clone();

  std::shared_ptr<double3DReg> vu =
      std::dynamic_pointer_cast<double3DReg>(storage);
  if (vu) return vu->clone();

  std::shared_ptr<double4DReg> vv =
      std::dynamic_pointer_cast<double4DReg>(storage);
  if (vv) return vv->clone();

  std::shared_ptr<double5DReg> vw =
      std::dynamic_pointer_cast<double5DReg>(storage);
  if (vw) return vw->clone();

  std::shared_ptr<double6DReg> vx =
      std::dynamic_pointer_cast<double6DReg>(storage);
  if (vx) return vx->clone();
  return nullptr;
}
