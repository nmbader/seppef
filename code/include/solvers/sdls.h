#ifndef sdls_H
#define sdls_H

#include "lsolver.h"

using namespace SEP;

// Linear steepest descent solver of min{norm2(Lm - d)}
class sdls : public lsolver{
protected:

public:

    // inherit all constructors from the base class
    using lsolver::lsolver;

    // destructor
    ~sdls(){}

    // run the solver for 1D problem
    void run(const oper1D * op, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat);

    // run the solver for 2D problem
    void run(const oper2D * op, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat);

    // test the solver for a full rank system Lm = d with 1D vectors
    void test();

};

#endif