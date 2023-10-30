#ifndef convolution2D_H
#define convolution2D_H

#include "oper2D.h"
#include "hypercube.h"
#include "float2DReg.h"

using namespace SEP;

// 2D convolution2D operator
class convolution2D : public oper2D{
protected:
    std::shared_ptr<float2DReg> _filter;
    bool _crossBounds;

public:
    // default constructor
    convolution2D(){_filter = nullptr;}

    // constructor
    convolution2D(const std::shared_ptr<float2DReg> filter, const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range){
        if(
            filter->getHyper()->getAxis(2).d != domain->getAxis(2).d ||
            filter->getHyper()->getAxis(1).d != domain->getAxis(1).d ||
            filter->getHyper()->getAxis(2).d != range->getAxis(2).d ||
            filter->getHyper()->getAxis(1).d != range->getAxis(1).d ||
            range->getAxis(2).o != filter->getHyper()->getAxis(2).o + domain->getAxis(2).o ||
            range->getAxis(1).o != filter->getHyper()->getAxis(1).o + domain->getAxis(1).o )     
                throw std::logic_error("Sampling mismatch or origins inconsistency for the convolution2D operator\n");

        _filter =  filter->clone();
        _domain = domain->clone();
        _range = range->clone();
        _crossBounds = true;
    }

    // destuctor
    ~convolution2D(){}

    // clone the object
    convolution2D * clone() const{
        convolution2D * conv = new convolution2D(_filter, _domain, _range);
        conv->crossBounds(_crossBounds);
        return conv;
    }

    // Define whether to roll the filter over the boundaries
    void crossBounds(const bool xBounds){
        _crossBounds = xBounds;
    }

    // forward operator: convolution
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator: correlation
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

#endif