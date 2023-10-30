#ifndef pef2D_H
#define pef2D_H

#include <stdexcept>
#include <cstdlib>
#include "float2DReg.h"
#include "cgls.h"

using namespace SEP;

//2D Prediction Error Filters
class pef2D {
protected:
    // leading 1 position on the first column
    // and size of gap following it if any
    unsigned int _lead, _gap;

public:
    // default constructor
    pef2D(){}

    // constructor
    pef2D(unsigned int lead=0, unsigned int gap=0){_lead = lead; _gap = gap;}

    // get the parameters
    unsigned int getLead(){return _lead;}

    unsigned int getGap(){return _gap;}

    // set the parameters
    void setLead(const unsigned int lead){_lead = lead;}

    void setGap(const unsigned int gap){_gap = gap;}

    // estimate stationary pef2D from data
    void estimate(std::shared_ptr<float2DReg> pef, const std::shared_ptr<float2DReg> dat, cgls &cg,
    const bool crossBounds = true, const unsigned int niter0 = 9999);

    // estimate stationary 2D PEF from data with mask
    void estimate(std::shared_ptr<float2DReg> pef, const std::shared_ptr<float2DReg> dat, 
                    const std::shared_ptr<float2DReg> mask, cgls &cg);

};

#endif