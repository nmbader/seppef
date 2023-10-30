#ifndef oper_H
#define oper_H

#include <stdexcept>
#include "hypercube.h"
#include "floatHyper.h"
#include "doubleHyper.h"
#include "float1DReg.h"
#include "double1DReg.h"
#include "float2DReg.h"
#include "double2DReg.h"
#include "float3DReg.h"
#include "double3DReg.h"
#include "float4DReg.h"
#include "double4DReg.h"

#ifdef DOUBLE_PRECISION
    typedef double data_t;
    typedef SEP::doubleHyper vecHyper;
    typedef SEP::double1DReg vec1DReg;
    typedef SEP::double2DReg vec2DReg;
    typedef SEP::double3DReg vec3DReg;
    typedef SEP::double4DReg vec4DReg;
#else
    typedef float data_t;
    typedef SEP::floatHyper vecHyper;
    typedef SEP::float1DReg vec1DReg;
    typedef SEP::float2DReg vec2DReg;
    typedef SEP::float3DReg vec3DReg;
    typedef SEP::float4DReg vec4DReg;
#endif

using namespace SEP;

// base class for linear operators
class oper {
protected:
    std::shared_ptr<hypercube> _domain;
    std::shared_ptr<hypercube> _range;

public:
    // default constructor
	oper(){}

    // virtual destructor
    virtual ~oper(){}

    // constructor
    oper(const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range){
        _domain = domain->clone();
        _range = range->clone();
    }

    // abstract method to clone the operator
    //virtual oper * clone() const = 0;

    // get domain
    std::shared_ptr<hypercube> getDomain() const {return _domain;}

    // get range
    std::shared_ptr<hypercube> getRange() const {return _range;}

    // set domain and range from hypercube
    void setDomainRange(const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range){
        _domain = domain->clone();
        _range = range->clone();
    }

    // set domain and range from floatHyper
    void setDomainRange(const std::shared_ptr<floatHyper> mod, const std::shared_ptr<floatHyper> dat){
        _domain = mod->getHyper()->clone();
        _range = dat->getHyper()->clone();
    }

    // check domain and range with hypercube
    inline void checkDomainRange(const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range) const {
        _domain->checkSame(domain);
        _range->checkSame(range);
    }

    // check domain and range with float2DReg
    inline void checkDomainRange(const std::shared_ptr<floatHyper> mod, const std::shared_ptr<floatHyper> dat) const {
        _domain->checkSame(mod->getHyper());
        _range->checkSame(dat->getHyper());
    }

};

#endif