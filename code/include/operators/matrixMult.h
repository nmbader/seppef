#ifndef matrixMult_H
#define matrixMult_H

#include "oper1D.h"
#include "float1DReg.h"
#include "float2DReg.h"
#include "floatHyperExt.h"
#include "hypercube.h"

using namespace SEP;

// Matrix multiplication with 1D vectors
class matrixMult : public oper1D{
private:
    std::shared_ptr<float2DReg> _mat;

public:
    // default constructor
    matrixMult(){_mat = nullptr;}

    // destructor
    ~matrixMult(){}

    // constructor from given size
    matrixMult(unsigned int nx, unsigned int nz){
        _mat = std::make_shared<float2DReg>(nz,nx);
        _domain = std::make_shared<hypercube>(axis(nx));
        _range = std::make_shared<hypercube>(axis(nz));
    }

    // constructor from float2DReg
    matrixMult(std::shared_ptr<float2DReg> mat){
        _mat = mat->clone();
        _domain = std::make_shared<hypercube>(mat->getHyper()->getAxis(2));
        _range = std::make_shared<hypercube>(mat->getHyper()->getAxis(1));
    }

    // cloning the object
    matrixMult * clone() const{
        matrixMult * mM = new matrixMult(_mat);
        return mM;
    }

    // fill the matrix with random numbers
    void random(float min =-1, float max=1){
        VECEXT::random(_mat, min, max);
    }

    // forward operator: matrix multiplication
    void forward(bool add, const std::shared_ptr<float1DReg> mod, std::shared_ptr<float1DReg> dat) const;

    // adjoint operator: matrix transpose multiplication
    void adjoint(bool add, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat) const;

};

#endif