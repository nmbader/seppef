#ifndef fkTransform_H
#define fkTransform_H

#include "complex2DReg.h"
#include "float2DReg.h"
#include "axis.h"
#include "hypercube.h"
#include <assert.h>


using namespace SEP;

// class to tranform from t-x (or z-x) to f-k space and vice versa
class fkTransform {
protected:
    std::shared_ptr<hypercube> _domain;
    std::shared_ptr<hypercube> _range;
    std::shared_ptr<hypercube> _crange;

public:
    // default constructor
	fkTransform(){}

    // destructor
    ~fkTransform(){}

    // constructor
    fkTransform(const std::shared_ptr<hypercube> domain){
        _domain = domain->clone();
        axis X = domain->getAxis(2);
        axis Z = domain->getAxis(1);
        Z.n = domain->getAxis(1).n/2 + 1;
        Z.d = 1.0/(domain->getAxis(1).n*domain->getAxis(1).d);
        Z.o = 0.0;
        X.d = 2*M_PI/(X.n*X.d);
        X.o = - X.d * floor((X.n-1)/2);
        _range = std::make_shared<hypercube>(Z, X);
        Z.n = domain->getAxis(1).n;
        Z.o = -(Z.n-1)/2 * Z.d;
        _crange = std::make_shared<hypercube>(Z, X);
    }

    // cloning the object
    fkTransform * clone() const{
        fkTransform * fk = new fkTransform(_domain);
        return fk;
    }

    // get domain
    std::shared_ptr<hypercube> getDomain() const {return _domain;}

    // get range for real domain
    std::shared_ptr<hypercube> getRange() const {return _range;}

    // get range for complex domain
    std::shared_ptr<hypercube> getCRange() const {return _crange;}

    // get number of frequencies
    unsigned int nf() const {return _range->getAxis(1).n;}

    // get frequency increment
    float df() const {return _range->getAxis(1).d;}

    // get spatial frequency increment
    float dk() const {return _range->getAxis(2).d;}

    // get origin of spatial frequency
    float ok() const {return _range->getAxis(2).o;}

    // set domain from hypercube and deduce range and crange
    void setDomainRange(const std::shared_ptr<hypercube> domain){
        _domain = domain->clone();
        axis X = domain->getAxis(2);
        axis Z = domain->getAxis(1);
        Z.n = domain->getAxis(1).n/2 + 1;
        Z.d = 1.0/(domain->getAxis(1).n*domain->getAxis(1).d);
        Z.o = 0.0;
        X.d = 2*M_PI/(X.n*X.d);
        X.o = - X.d * floor((X.n-1)/2);
        _range = std::make_shared<hypercube>(Z, X);
        Z.n = domain->getAxis(1).n;
        Z.o = -(Z.n-1)/2 * Z.d;
        _crange = std::make_shared<hypercube>(Z, X);
    }

    // set domain from float2DReg and deduce range and crange
    void setDomainRange(const std::shared_ptr<float2DReg> mod){
        _domain = mod->getHyper()->clone();
        axis X = _domain->getAxis(2);
        axis Z = _domain->getAxis(1);
        Z.n = _domain->getAxis(1).n/2 + 1;
        Z.d = 1.0/(_domain->getAxis(1).n*_domain->getAxis(1).d);
        Z.o = 0.0;
        X.d = 2*M_PI/(X.n*X.d);
        X.o = - X.d * floor((X.n-1)/2);
        _range = std::make_shared<hypercube>(Z, X);
        Z.n = _domain->getAxis(1).n;
        Z.o = -(Z.n-1)/2 * Z.d;
        _crange = std::make_shared<hypercube>(Z, X);
    }

    // check domain and range for real domain
    inline void checkDomainRange(const std::shared_ptr<float2DReg> mod, const std::shared_ptr<complex2DReg> dat) const {
        _domain->checkSame(mod->getHyper());
        _range->checkSame(dat->getHyper());
    }

    // check domain and range for complex domain
    inline void checkDomainRange(const std::shared_ptr<complex2DReg> mod, const std::shared_ptr<complex2DReg> dat) const {
        _domain->checkSame(mod->getHyper());
        _crange->checkSame(dat->getHyper());
    }

    // forward fkTransform
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<complex2DReg> dat) const;

    // inverse fkTransform
    void inverse(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<complex2DReg> dat) const;

    // forward fkTransform on complex data
    void cforward(bool add, const std::shared_ptr<complex2DReg> mod, std::shared_ptr<complex2DReg> dat) const;

    // inverse fkTransform on complex data
    void cinverse(bool add, std::shared_ptr<complex2DReg> mod, const std::shared_ptr<complex2DReg> dat) const;

    // test the forward-inverse
    void testInverse();

    // test the cforward-cinverse
    void testCInverse();

};

#endif