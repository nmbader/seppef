#ifndef oper1D_H
#define oper1D_H

#include <stdexcept>
#include "hypercube.h"
#include "float1DReg.h"
#include "oper.h"

using namespace SEP;

// base class for 1D linear operators
class oper1D : public oper {

public:
    // default constructor
	oper1D(){}

    // virtual destructor
    virtual ~oper1D(){}

    // abstract method to clone the operator
    virtual oper1D * clone() const = 0;

    // abstract forward operator
    virtual void forward(bool add, const std::shared_ptr<float1DReg> mod, std::shared_ptr<float1DReg> dat) const = 0;

    // abstract adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat) const = 0;

    // perform the dot product test of the operator
    void dotProduct();

};

#endif