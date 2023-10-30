#ifndef cgls_H
#define cgls_H

#include "lsolver.h"

using namespace SEP;

// Linear conjugate gradient least squared CGLS solver of min{norm2(Lm - d)}
// according to Aster (Parameter Estimation and Inverse Problems), chap 6
class cgls : public lsolver{
protected:

public:

    // inherit all constructors from the base class
    using lsolver::lsolver;

    // destructor
    ~cgls(){}

    // run the solver for 1D problem
    void run(const oper1D * op, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat);

    // run the solver for 2D problem
    void run(const oper2D * op, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat);

    
};

#endif