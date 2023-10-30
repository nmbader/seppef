#include <fstream>
#include "lsolverReg.h"
#include "matrixMult.h"
#include "identity1D.h"

using namespace SEP;

void lsolverReg::test(){

    matrixMult mult(50,200);
    std::shared_ptr<float1DReg> mod (new float1DReg(50));
    std::shared_ptr<float1DReg> mod0 (new float1DReg(50));
    std::shared_ptr<float1DReg> mod1 (new float1DReg(50));
    std::shared_ptr<float1DReg> dat (new float1DReg(200));

    mult.random(-1,1);
    VECEXT::random(mod,-10,10);
    VECEXT::random(mod0,-10,10);
    VECEXT::random(mod1,-10,10);
    mult.forward(false,mod,dat);

    identity1D Id (mod->getHyper());
    float eps0 = 1e-03;
    float eps1 = 10;

    this->run(&mult, &Id, eps0, mod0, dat);
    this->run(&mult, &Id, eps1, mod1, dat);

    // L2 squared error of estimated model
    double trueNorm, estimatedNorm0, estimatedNorm1, dataNorm, error0, error1;
    dataNorm = VECEXT::norm2(dat);
    trueNorm = VECEXT::norm2(mod);
    estimatedNorm0 = VECEXT::norm2(mod0);
    estimatedNorm1 = VECEXT::norm2(mod1);
    mod->scale(-1);
    mod0->add(mod);
    mod1->add(mod);
    error0 = VECEXT::norm2(mod0);
    error1 = VECEXT::norm2(mod1);

    std::clog << "data norm squared = " << dataNorm << std::endl;
    std::clog << "true model norm squared = " << trueNorm << std::endl;
    std::clog << "estimated model norm squared (weak reg) = " << estimatedNorm0 << std::endl;
    std::clog << "estimated model norm squared (strong reg) = " << estimatedNorm1 << std::endl;
    std::clog << "model error squared (weak reg) = " << error0 << std::endl;
    std::clog << "model error squared (strong reg) = " << error1 << std::endl;
}