#ifndef cglsSuper_H
#define cglsSuper_H

#include "float1DReg.h"
#include "float2DReg.h"
#include "oper.h"
#include "oper2D.h"
#include "oper1D.h"
#include "cgls.h"


using namespace SEP;

// CGLS solver of min{norm2(L1m - d1) + norm2(L2m - d2)}
class cglsSuper : public cgls{

public:
    // Inheriting all constructors from CGLS
    using cgls::cgls;

    // run the solver for 1D problem
    void run(const oper1D * L1, const oper1D * L2, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat1, const std::shared_ptr<float1DReg> dat2);

    // run the solver for 2D problem
    void run(const oper2D * L1, const oper2D * L2, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat1, const std::shared_ptr<float2DReg> dat2);

    // test the solver for a full rank system L1m = d1, L2m = d2 with 1D vectors
    void test();
    
};

#endif