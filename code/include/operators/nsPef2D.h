#ifndef nsPef2D_H
#define nsPef2D_H

#include "oper2D.h"
#include "hypercube.h"
#include "float2DReg.h"
#include "floatHyperExt.h"
#include "cgls.h"

using namespace SEP;

// Non-stationary 2D PEF operator
class nsPef2D : public oper2D{
protected:
    // The initial PEF; serves also for the structure of the filter
    std::shared_ptr<float2DReg> _pef;

    // leading 1 position on the first column
    unsigned int _lead;

public:
    // default constructor
    nsPef2D(){}

    // constructor
    nsPef2D(const std::shared_ptr<float2DReg> pef,
    const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range,
    unsigned int lead=0){
  
        // check if domain and range are the same
        domain->checkSame(range);

        if (lead > pef->getHyper()->getAxis(1).n-1)
            throw std::logic_error("The lead index exceeds nz");

        if ((pef->getHyper()->getAxis(2).n > domain->getAxis(2).n) || (pef->getHyper()->getAxis(1).n > domain->getAxis(1).n))
            throw std::logic_error("data is smaller than filter");
        

std::clog << "The PEF will be reset to be compliant with the domain and the leading 1\n";
        
        // set the hypercube for the PEF, set originx to zero and originz according to the leading 1
        _pef =  pef->clone();
        axis X = domain->getAxis(2);
        axis Z = domain->getAxis(1);
        X.o = 0;
        Z.o = -Z.d * lead;
        std::shared_ptr<hypercube> hyper (new hypercube(Z, X));
        VECEXT::copyHyper(_pef, hyper);

        // set the pef to zero before the lead and to 1 at the lead
        for (int i=0; i<lead; i++)
            (*_pef->_mat)[0][i] = 0;
        (*_pef->_mat)[0][lead] = 1;

        
        _domain = domain->clone();
        _range = range->clone();
        _lead = lead;

    }

    // destuctor
    virtual ~nsPef2D(){}

    // get the parameters
    unsigned int getLead(){return _lead;}

    // set the parameters
    void setPef(const std::shared_ptr<float2DReg> pef){_pef = pef;}
    void setLead(const unsigned int lead){_lead = lead;}


    // forward and adjoint operators are based on the version selected. All versions loop over the time axis first,
    // and away from boundaries

    // forward operator: convolution
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const=0;

    // adjoint operator: correlation
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const=0;

};


// Jon Claerbout's style with fixed step length
// user defined or calculated using 'computeDatNorm2'
class nsPef2D_v0 : public nsPef2D{
protected:
    // The data from which to estimate the PEF
    std::shared_ptr<float2DReg> _dat0;

    // step length for v0 class
    // it is used as damping factor for step length calculation for v1 class
    double _epsilon;

    // exposed forwarding constructor
    nsPef2D_v0(const std::shared_ptr<float2DReg> pef,
    const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range,
    unsigned int lead = 0):nsPef2D(pef, domain, range, lead){}

public:
    // Final PEF after running forward or adjoint operator
    mutable std::shared_ptr<float2DReg> _peff;

    // default constructor
    nsPef2D_v0(){}

    // constructor
    nsPef2D_v0(const std::shared_ptr<float2DReg> pef, const std::shared_ptr<float2DReg> dat0,
    const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range,
    unsigned int lead=0, double epsilon = 0):nsPef2D(pef, domain, range, lead){
  
        // check if domain and range are the same as the data from which the PEF is estimated
        domain->checkSame(dat0->getHyper());

        _dat0 = dat0->clone();
        
        _epsilon = epsilon;

        if (_epsilon == 0){
std::clog << "Epsilon value will be set based on data norm2 squared and PEF structure = ";
            std::shared_ptr<float2DReg> Norm2 = this->computeDatNorm2();
            this->setEpsilon(1/Norm2->max());
std::clog << _epsilon << "\n";
        }
    }

    // destuctor
    ~nsPef2D_v0(){}

    // clone the object
    nsPef2D_v0 * clone() const{
        nsPef2D_v0 * nsPef = new nsPef2D_v0(_pef, _dat0, _domain, _range, _lead, _epsilon);
        return nsPef;
    }

    // set the parameters
    void setEpsilon(const double epsilon){_epsilon = epsilon;}

    // get the parameters
    double getEpsilon(){return _epsilon;}

    // compute the data L2 norm squared from which the PEF is estimated
    std::shared_ptr<float2DReg> computeDatNorm2();

    // forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};



// Jon Claerbout's style with variable step length
// calculated automatically and damped by epsilon: SL = r2/(r2.d2 + eps)
class nsPef2D_v1 : public nsPef2D_v0{
protected:

public:
    // inherit all constructors from the base class
    using nsPef2D_v0::nsPef2D_v0;

    // destuctor
    ~nsPef2D_v1(){}

    // clone the object
    nsPef2D_v1 * clone() const{
        nsPef2D_v1 * nsPef = new nsPef2D_v1(_pef, _dat0, _domain, _range, _lead, _epsilon);
        return nsPef;
    }

    // forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};



// Jon Claerbout's style with variable step length
// calculated automatically and damped by epsilon: SL = r2/(r2.d2 + eps)
// and with weighting in x and z
class nsPef2D_v2 : public nsPef2D_v0{
protected:
    float _weightx, _weightz;

public:
    // inherit all constructors from the base class
    using nsPef2D_v0::nsPef2D_v0;

    // destuctor
    ~nsPef2D_v2(){}

    // clone the object
    nsPef2D_v2 * clone() const{
        nsPef2D_v2 * nsPef = new nsPef2D_v2(_pef, _dat0, _domain, _range, _lead, _epsilon);
        nsPef->setWeightx(_weightx);
        return nsPef;
    }

    // set the weights
    void setWeightx(const float wx = 0){
        
        assert(wx >= 0);
        assert(wx <= 1);

        _weightx = wx;
        _weightz = 1 - wx;
    }

    // forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

// Jon Claerbout's style with variable step length
// calculated automatically: SL = 1/(eps.d2 + (1-eps).d2_max)
// and with weighting in x and z
class nsPef2D_v3 : public nsPef2D_v0{
protected:
    float _weightx, _weightz;
    double _d2max;

public:
    // default constructor
    nsPef2D_v3(){}

    // constructor
    nsPef2D_v3(const std::shared_ptr<float2DReg> pef, const std::shared_ptr<float2DReg> dat0,
    const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range,
    unsigned int lead=0, float epsilon = 0, float wx = 0):nsPef2D_v0(pef, domain, range, lead){
  
        // check if domain and range are the same as the data from which the PEF is estimated
        domain->checkSame(dat0->getHyper());

        _dat0 = dat0->clone();
        
        assert((epsilon >= 0) && (epsilon <= 1));
        assert((wx >= 0) && (wx <= 1));

        _epsilon = epsilon;
        _weightx = wx;
        _weightz = 1 - wx;

        std::shared_ptr<float2DReg> Norm2 = this->computeDatNorm2();
        _d2max = Norm2->max();
    }

    // destuctor
    ~nsPef2D_v3(){}

    // clone the object
    nsPef2D_v3 * clone() const{
        nsPef2D_v3 * nsPef = new nsPef2D_v3(_pef, _dat0, _domain, _range, _lead, _epsilon);
        nsPef->setWeightx(_weightx);
        return nsPef;
    }

    // set the weights
    void setWeightx(const float wx = 0){
        
        assert(wx >= 0);
        assert(wx <= 1);

        _weightx = wx;
        _weightz = 1 - wx;
    }

    // forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};


// Sergey Fomel's + Jon Claerbout's style with variable step length
// calculated automatically and damped by epsilon: SL = 1.0/(d2 + eps)
class nsPef2D_v4 : public nsPef2D_v0{
protected:

public:
    // inherit all constructors from the base class
    using nsPef2D_v0::nsPef2D_v0;

    // destuctor
    ~nsPef2D_v4(){}

    // clone the object
    nsPef2D_v4 * clone() const{
        nsPef2D_v4 * nsPef = new nsPef2D_v4(_pef, _dat0, _domain, _range, _lead, _epsilon);
        return nsPef;
    }

    // forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};

// Jon Claerbout's style with variable step length
// calculated automatically and damped by epsilon: SL = 1/(d2 + eps/r2)
// and with weighting in x and z
class nsPef2D_v5 : public nsPef2D_v0{
protected:
    float _weightx, _weightz;

public:
    // inherit all constructors from the base class
    using nsPef2D_v0::nsPef2D_v0;

    // destuctor
    ~nsPef2D_v5(){}

    // clone the object
    nsPef2D_v5 * clone() const{
        nsPef2D_v5 * nsPef = new nsPef2D_v5(_pef, _dat0, _domain, _range, _lead, _epsilon);
        nsPef->setWeightx(_weightx);
        return nsPef;
    }

    // set the weights
    void setWeightx(const float wx = 0){
        
        assert(wx >= 0);
        assert(wx <= 1);

        _weightx = wx;
        _weightz = 1 - wx;
    }

    // forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};


// Non-stat PEF by patch where stationary PEF is estimated
// in each local window by steepest descent
// PEF interpolation is possible between windows
class nsPef2D_v9 : public nsPef2D{
protected:
    // The data from which to estimate the PEF
    std::shared_ptr<float2DReg> _dat0;

    // windows size for the local patches
    int _winx, _winz;

    // number of iterations to be used for steepest descent for PEF estimation
    unsigned int _niter;

public:
    // Final PEF after running forward or adjoint operator
    mutable std::shared_ptr<float2DReg> _peff;

    // default constructor
    nsPef2D_v9(){}

    // constructor
    nsPef2D_v9(const std::shared_ptr<float2DReg> pef, const std::shared_ptr<float2DReg> dat0,
    const std::shared_ptr<hypercube> domain, const std::shared_ptr<hypercube> range,
    unsigned int lead=0, int winx=10, int winz=10,
    unsigned int niter=5):nsPef2D(pef, domain, range, lead){
  
        // check if domain and range are the same as the data from which the PEF is estimated
        domain->checkSame(dat0->getHyper());

        _dat0 = dat0->clone();
        
        _winx = winx;
        _winz = winz;
        _niter = niter;

    }

    // destuctor
    ~nsPef2D_v9(){}

    // clone the object
    nsPef2D_v9 * clone() const{
        nsPef2D_v9 * nsPef = new nsPef2D_v9(_pef, _dat0, _domain, _range, _lead, _winx, _winz, _niter);
        return nsPef;
    }

    // forward operator
    virtual void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    virtual void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

};


#endif