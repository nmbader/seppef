#ifndef zeroPhaseFiltering2D_H
#define zeroPhaseFiltering2D_H

#include "oper2D.h"
#include "float1DReg.h"
#include "float2DReg.h"

using namespace SEP;

// zero phase filtering on 2D data (1D convolution)
class zeroPhaseFiltering2D : public oper2D {
     // Filter must be 1D and have odd number of samples
     std::shared_ptr<float1DReg> _filter;

public:
    // default constructor
    zeroPhaseFiltering2D(){_filter = nullptr;}

    // destructor
    ~zeroPhaseFiltering2D(){}

    // constructor
    zeroPhaseFiltering2D(const std::shared_ptr<float1DReg> filter, const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range) {

        if(2*(filter->getHyper()->getAxis(1).n/2) == filter->getHyper()->getAxis(1).n)     
            throw std::logic_error("Number of filter samples is not odd. Cannot be used for symmetric convolution.\n");

        domain->checkSame(range);

        if(
            filter->getHyper()->getAxis(1).d != domain->getAxis(1).d ||
            std::abs(range->getAxis(1).o - filter->getHyper()->getAxis(1).o - filter->getHyper()->getAxis(1).d*(int)(filter->getHyper()->getAxis(1).n/2)) > 1e-06 )     
                throw std::logic_error("Sampling mismatch or origins inconsistency");

        if (filter->getHyper()->getAxis(1).n > range->getAxis(1).n)
            throw std::logic_error("filter longer than data\n");

        _filter = filter->clone();
        _domain = domain->clone();
        _range = range->clone();
    }

    // cloning the object
    zeroPhaseFiltering2D * clone() const {
        zeroPhaseFiltering2D * conv = new zeroPhaseFiltering2D(_filter, _domain, _range);
        return conv;
    }

    // forward operator: convolution
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator: correlation
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

#endif