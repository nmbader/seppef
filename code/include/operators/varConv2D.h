#ifndef varConv2D_H
#define varConv2D_H

#include <assert.h>

#include "oper2D.h"
#include "hypercube.h"
#include "float2DReg.h"

using namespace SEP;

// 2D convolution operator with variable filters defined on a grid
class varConv2D : public oper2D{
protected:
    // size of the local filter
    int _nx, _nz;

    // leading 1 position on the first column
    unsigned int _lead;

    // grid of filters
    int _xinc, _zinc; 
    int _nfilx, _nfilz;

    // all filters
    std::shared_ptr<float2DReg> _allfil;

public:
    // default constructor
    varConv2D(){}

    // constructor
    varConv2D(const std::shared_ptr<float2DReg> filters, const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range,
    int nx, int nz, int xinc, int zinc, unsigned int lead=0){
        
        domain->checkSame(range);

        _domain = domain->clone();
        _range = range->clone();

        _nx = nx;
        _nz = nz;
        _xinc = xinc;
        _zinc = zinc;
        _lead = lead;
        _nfilx = ceil(domain->getAxis(2).n / (float)xinc) + 1;
        _nfilz = ceil(domain->getAxis(1).n / (float)zinc) + 1;

        assert(_nfilz * _nz == filters->getHyper()->getAxis(1).n);
        assert(_nfilx * _nx == filters->getHyper()->getAxis(2).n);

        _allfil = filters->clone();

    }

    // destuctor
    ~varConv2D(){}

    // clone the object
    varConv2D * clone() const{
        varConv2D * conv = new varConv2D(_allfil, _domain, _range, _nx, _nz, _xinc, _zinc, _lead);
        return conv;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

#endif
