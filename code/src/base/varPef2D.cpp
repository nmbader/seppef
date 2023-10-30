#include "varPef2D.h"
#include "varDatConv2D.h"
#include "selection2D.h"
#include "chainOper2D.h"

void varPef2D::estimate(const std::shared_ptr<float2DReg> dat, cgls &cg){

    _domain->checkSame(dat->getHyper());

    // Solving D.Ku.p = -D.Kn.p
    // D: data convolution operator
    // Kn, Ku: selector operators
    // p: variable PEF to estimate (= _allpef)

    // build the data convolution operator with variable filter
    axis X, Z;
    X.n = _nx;
    Z.n = _nz;
    Z.d = _domain->getAxis(1).d;
    X.d = _domain->getAxis(2).d;
    X.o = 0;
    Z.o = -Z.d * _lead;
    std::shared_ptr<hypercube> hyper (new hypercube(Z, X));
    varDatConv2D D(dat, hyper, _xinc, _zinc, _lead);

    // build the selection operators to separate PEF known and unknown coefficients
    std::shared_ptr<float2DReg> selector = _allpef->clone();
    selector->zero();

    for (int ix=0; ix<_nfilx; ix++){
        for (int iz=0; iz<_nfilz; iz++){
            for (int izf=0; izf <= _lead; izf++){
                (*selector->_mat)[ix*_nx][iz*_nz+izf] = 1;
            }
        }
    }

    selection2D Kn(selector);
    selection2D Ku(selector);
    Ku.reverse();

    // build the residual vector res = -D.Kn.p
    std::shared_ptr<float2DReg> res (new float2DReg(dat->getHyper()));
    std::shared_ptr<float2DReg> KnP (new float2DReg(_allpef->getHyper()));
    Kn.forward(false, _allpef, KnP);

    D.forward(false, KnP, res);
    res->scale(-1);

    // build the chain operator D.Ku
    chainOper2D DKu(&D, &Ku);

    // initialize and run the CGLS solver
    unsigned int niter = (_nx*_nz - _lead-1)*_nfilx*_nfilz;
    niter = std::min(niter, cg.getNiter());
    cg.setNiter(niter);

std::clog << "Maximum number of iterations for CG solver set to: " <<niter<< "\n";

    cg.run(&DKu, _allpef, res);

}