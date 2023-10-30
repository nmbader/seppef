#include <fstream>
#include "lsolver.h"
#include "matrixMult.h"
#include "floatHyperExt.h"

using namespace SEP;

void lsolver::test(){

    matrixMult mult(50,200);
    std::shared_ptr<float1DReg> mod (new float1DReg(50));
    std::shared_ptr<float1DReg> mod0 (new float1DReg(50));
    std::shared_ptr<float1DReg> dat (new float1DReg(200));

    mult.random(-1,1);
    VECEXT::random(mod,-10,10);
    mult.forward(false,mod,dat);

    this->run(&mult, mod0, dat);

    // L2 squared error of estimated model
    double trueNorm, estimatedNorm, dataNorm, error;
    trueNorm = VECEXT::norm2(mod);
    estimatedNorm = VECEXT::norm2(mod0);
    dataNorm = VECEXT::norm2(dat);
    mod0->scale(-1);
    mod->add(mod0);
    error = VECEXT::norm2(mod);

    std::clog << "data norm squared = " << dataNorm << std::endl;
    std::clog << "true model norm squared = " << trueNorm << std::endl;
    std::clog << "estimated model norm squared = " << estimatedNorm << std::endl;
    std::clog << "model error squared = " << error << std::endl;
}