#ifndef expVarConv2D_H
#define expVarConv2D_H

#include "varConv2D.h"

using namespace SEP;

// 2D convolution operator with variable filters defined on a grid, expanded when applied
class expVarConv2D : public varConv2D{

public:
    // inherit all constructors from the base class
    using varConv2D::varConv2D;

    // destuctor
    ~expVarConv2D(){}

    // forward and adjoint operators are the same as for the base class, except that the filters
    // are expanded before application

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

#endif
