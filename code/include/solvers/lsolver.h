#ifndef lsolver_H
#define lsolver_H

#include "float1DReg.h"
#include "float2DReg.h"
#include "oper.h"
#include "oper1D.h"
#include "oper2D.h"

using namespace SEP;

// Iterative solver of min{norm2(Lm - d)}
class lsolver{
protected:
    unsigned int _niter;
    float _threshold; // stopper based on rate of convergence

public:
    std::vector <float> _func;
    std::shared_ptr<floatHyper> _grad;
    std::shared_ptr<floatHyper> _res;

    // default constructor
    lsolver(){}

    // constructor
    lsolver(unsigned int niter, float threshold=0){
        _niter = niter;
        _threshold = threshold;
        _func = {};
        _grad = nullptr;
        _res = nullptr;
    }

    // destructor
    ~lsolver(){}

    // get the parameters
    unsigned int getNiter() const {return _niter;}
    float getThreshold() const {return _threshold;}

    // set the number of iterations
    void setNiter(unsigned int niter) {_niter = niter;}

    // set the threshold
    void setThreshold(float threshold){_threshold = threshold;}

    // run the solver for 1D problem
    virtual void run(const oper1D * op, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat) = 0;

    // run the solver for 2D problem
    virtual void run(const oper2D * op, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) = 0;

    // reset the objective function
    void resetFunc(){
        _func.clear();
    }

    // test the solver for a full rank system Lm = d with 1D vectors
    void test();
};

#endif