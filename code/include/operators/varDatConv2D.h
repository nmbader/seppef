#ifndef varDatConv2D_H
#define varDatConv2D_H

#include "axis.h"
#include "oper2D.h"
#include "hypercube.h"
#include "float2DReg.h"
#include "floatHyperExt.h"

using namespace SEP;

// Operator to convolve data with variable filter defined on a grid: D.f
class varDatConv2D : public oper2D{
protected:
    std::shared_ptr<float2DReg> _dat0;
    std::shared_ptr<hypercube> _fil0;
    int _xinc, _zinc;
    int _nfilx, _nfilz;
    int _lead;

public:
    // default constructor
    varDatConv2D(){_dat0 = nullptr; _fil0 = nullptr;}

    // constructor
    varDatConv2D(const std::shared_ptr<float2DReg> dat0, const std::shared_ptr<hypercube> fil0,
    int xinc, int zinc, int lead=0){

        assert(lead < fil0->getAxis(1).n);

        if ((fil0->getAxis(2).n > dat0->getHyper()->getAxis(2).n) || (fil0->getAxis(1).n > dat0->getHyper()->getAxis(1).n))
            throw std::logic_error("data is smaller than filter");

        _range = dat0->getHyper()->clone();

        _xinc = xinc;
        _zinc = zinc;
        _lead = lead;
        _nfilx = ceil(_range->getAxis(2).n / (float)xinc) + 1;
        _nfilz = ceil(_range->getAxis(1).n / (float)zinc) + 1;


std::clog << "The filter will be reset to be compliant with the domain and the leading sample\n";
        
        // set the hypercube for the filter, set originx to zero and originz according to the leading sample
        axis Z = fil0->getAxis(1);
        axis X = fil0->getAxis(2);
        Z.d = _range->getAxis(1).d;
        X.d = _range->getAxis(2).d;
        X.o = 0;
        Z.o = -Z.d * lead;
        std::shared_ptr<hypercube> hyper (new hypercube(Z, X));
        _fil0 = hyper;

        Z.o = 0;
        Z.n = _nfilz * Z.n;
        X.n = _nfilx * X.n;

        _domain = std::make_shared<hypercube>(Z,X);
        _dat0 = dat0->clone();
    }

    // destuctor
    ~varDatConv2D(){}

    // clone the object
    varDatConv2D * clone() const{
        varDatConv2D * conv = new varDatConv2D(_dat0, _fil0, _xinc, _zinc, _lead);
        return conv;
    }

    // get the parameters
    int getNfilx() const {return _nfilx;}
    int getNfilz() const {return _nfilz;}

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

#endif