#ifndef cglsReg_H
#define cglsReg_H

#include "lsolverReg.h"

using namespace SEP;

// CGLS solver of min{norm2(Lm - d) + eps2.norm2(D(m - m0))}
class cglsReg : public lsolverReg{

public:
    // Inheriting all constructors from base class
    using lsolverReg::lsolverReg;

    // destructor
    ~cglsReg(){}

    // run the solver for 1D problem
    void run(const oper1D * L, const oper1D * D, const float eps, std::shared_ptr<float1DReg> mod,
                    const std::shared_ptr<float1DReg> dat, std::shared_ptr<float1DReg> mod0 = nullptr);

    // run the solver for 2D problem
    void run(const oper2D * L, const oper2D * D, const float eps, std::shared_ptr<float2DReg> mod,
                    const std::shared_ptr<float2DReg> dat, std::shared_ptr<float2DReg> mod0 = nullptr);
    
};

#endif