#ifndef varPef2D_H
#define varPef2D_H

#include <assert.h>
#include "float2DReg.h"
#include "cgls.h"

using namespace SEP;

// Variable 2D Prediction Error Filters
// Bob Clapp's style with multiple PEFs defined on a regular grid
// The application of the PEF is by 4 pts linear interpolation
class varPef2D {
protected:
    // size of the local PEF
    int _nx, _nz;

    // leading 1 position on the first column
    unsigned int _lead;

    // grid of PEFs
    int _xinc, _zinc; 
    int _nfilx, _nfilz;

    // data domain form which to estimate the filter
    std::shared_ptr<hypercube> _domain;

public:
    // list of all PEFs defined on the grid, initialized to the same PEF provided
    std::shared_ptr<float2DReg> _allpef;

    // default constructor
    varPef2D(){}

    // constructor
    varPef2D(const std::shared_ptr<hypercube> domain, int nx, int nz, 
    int xinc, int zinc, int lead=0, std::shared_ptr<float2DReg> allpef=nullptr){

        assert(lead < nz);

        if ((nx > domain->getAxis(2).n) || (nz > domain->getAxis(1).n))
            throw std::logic_error("data is smaller than filter");

        _nx = nx;
        _nz = nz;
        _xinc = xinc;
        _zinc = zinc;
        _lead = lead;
        _nfilx = ceil(domain->getAxis(2).n / (float)xinc) + 1;
        _nfilz = ceil(domain->getAxis(1).n / (float)zinc) + 1;
        _domain = domain->clone();
        
        // initialize the filters _allpef
        axis Z = domain->getAxis(1);
        axis X = domain->getAxis(2);
        X.o = 0;
        Z.o = 0;
        Z.n = _nfilz * nz;
        X.n = _nfilx * nx;

        
        if (allpef != nullptr){
            std::shared_ptr<hypercube> hyper (new hypercube(Z,X));
            allpef->getHyper()->checkSame(hyper);
            _allpef = allpef->clone();
        }
        else{
            _allpef = std::make_shared<float2DReg>(Z,X);
            _allpef->zero();
        }

        for (int ix=0; ix<_nfilx; ix++){
           for (int iz=0; iz<_nfilz; iz++){
               for (int i=0; i<lead; i++){
                   (*_allpef->_mat)[ix*nx][iz*nz+i] = 0;
               }
                (*_allpef->_mat)[ix*nx][iz*nz+lead] = 1;
           }
        }
    }

    // destuctor
    ~varPef2D(){}

    // get the parameters
    unsigned int getLead(){return _lead;}

    // populate _allpef with copies of a given PEF
    void populate(const std::shared_ptr<float2DReg> pef){
        
        assert(_nz == pef->getHyper()->getAxis(1).n);
        assert(_nx == pef->getHyper()->getAxis(2).n);

        std::shared_ptr<float2DReg> pef0 = pef->clone();

        for (int iz=0; iz<_lead; iz++){
            (*pef0->_mat)[0][iz] = 0;
        }
        (*pef0->_mat)[0][_lead] = 1;

        for (int ix=0; ix<_nfilx; ix++){
           for (int iz=0; iz<_nfilz; iz++){
               for (int ixf=0; ixf<_nx; ixf++){
                   for (int izf=0; izf<_nz; izf++){
                       (*_allpef->_mat)[ix*_nx+ixf][iz*_nz+izf] = (*pef0->_mat)[ixf][izf];
                   }
               }
           }
       }
    }

    // reset the grid of PEFs with only ones at the leading index
    void reset(){

        _allpef->zero();

        for (int ix=0; ix<_nfilx; ix++){
           for (int iz=0; iz<_nfilz; iz++){
                (*_allpef->_mat)[ix*_nx][iz*_nz+_lead] = 1;
           }
        }
    }

    // estimate the _allpef from data
    void estimate(const std::shared_ptr<float2DReg> dat, cgls &cg);

};

#endif