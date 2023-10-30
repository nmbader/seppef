#ifndef oper2D_H
#define oper2D_H

#include <stdexcept>
#include "hypercube.h"
#include "float2DReg.h"
#include "oper.h"


using namespace SEP;

// base class for 2D linear oper2Dators
class oper2D : public oper {

public:
    // default constructor
	oper2D(){}

    // virtual destructor
    virtual ~oper2D(){}

    // abstract method to clone the operator
    virtual oper2D * clone() const = 0;

    // abstract forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const = 0;

    // abstract adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const = 0;

    // perform the dot product test of the operator
    void dotProduct();


};

#endif