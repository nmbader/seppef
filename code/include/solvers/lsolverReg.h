#ifndef lsolverReg_H
#define lsolverReg_H

#include "float1DReg.h"
#include "float2DReg.h"
#include "oper.h"
#include "oper1D.h"
#include "oper2D.h"

using namespace SEP;

// Iterative solver of min{norm2(Lm - d) + eps2.norm2(D(m - m0))}
class lsolverReg{
protected:
    unsigned int _niter;
    float _threshold; // stopper based on rate of convergence

public:
    std::vector <float> _func;
    std::shared_ptr<floatHyper> _grad;
    std::shared_ptr<floatHyper> _m_grad;
    std::shared_ptr<floatHyper> _d_res;
    std::shared_ptr<floatHyper> _m_res;

    // default constructor
    lsolverReg(){}

    // constructor
    lsolverReg(unsigned int niter, float threshold=0){
        _niter = niter;
        _threshold = threshold;
        _func = {};
        _grad = nullptr;
        _m_grad = nullptr;
        _d_res = nullptr;
        _m_res = nullptr;
    }

    // destructor
    ~lsolverReg(){}

    // get the parameters
    unsigned int getNiter() const {return _niter;}
    float getThreshold() const {return _threshold;}

    // set the parameters
    void setNiter(unsigned int niter) {_niter = niter;}
    void setThreshold(float threshold){_threshold = threshold;}

    // run the solver for 1D problem
    virtual void run(const oper1D * L, const oper1D * D, const float eps, std::shared_ptr<float1DReg> mod,
                    const std::shared_ptr<float1DReg> dat, std::shared_ptr<float1DReg> mod0 = nullptr) = 0;

    // run the solver for 2D problem
    virtual void run(const oper2D * L, const oper2D * D, const float eps, std::shared_ptr<float2DReg> mod,
                    const std::shared_ptr<float2DReg> dat, std::shared_ptr<float2DReg> mod0 = nullptr) = 0;

    // reset the objective function
    void resetFunc(){
        _func.clear();
    }

    // test the solver for a full rank system Lm = d, eps << or >>, D=Id, mod0 = 0 with 1D vectors
    void test();
    
};

#endif