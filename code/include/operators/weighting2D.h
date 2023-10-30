#ifndef weighting2D_H
#define weighting2D_H

#include "oper2D.h"

using namespace SEP;

// weighting operator (domain = range)
class weighting2D : public oper2D{
private:
    std::shared_ptr<float2DReg> _weights;

public:
    // default constructor
    weighting2D(){_weights = nullptr;}

    // destructor
    ~weighting2D(){}

    // constructor
    weighting2D(const std::shared_ptr<float2DReg> weights){
        _weights = weights->clone();
        _domain = weights->getHyper()->clone();
        _range = weights->getHyper()->clone();
    }

    // cloning the object
    weighting2D * clone() const{
        weighting2D * wei = new weighting2D(_weights);
        return wei;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

    // inverse the weights
    void inverse(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

#endif