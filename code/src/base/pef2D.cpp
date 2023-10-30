#include "pef2D.h"
#include "hypercube.h"
#include "axis.h"
#include "floatHyperExt.h"
#include "selection2D.h"
#include "convolution2D.h"
#include "chainOper2D.h"


void pef2D::estimate(std::shared_ptr<float2DReg> pef, const std::shared_ptr<float2DReg> dat, cgls &cg,
const bool crossBounds, const unsigned int niter0){

    if (_lead > pef->getHyper()->getAxis(1).n-1)
            throw std::logic_error("The lead index exceeds nz");

    if (_lead + _gap > pef->getHyper()->getAxis(1).n-1) 
        throw std::logic_error("Gap too large or filter too short");

    if ((VECEXT::getnx(pef) > VECEXT::getnx(dat)) || (VECEXT::getnz(pef) > VECEXT::getnz(dat)))
        throw std::logic_error("data is smaller than filter");

    // copy the hypercube from data to the PEF without the size and set originx to zero and originz according to the leading 1
    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);
    X.o = 0;
    Z.o = -Z.d * _lead;
    std::shared_ptr<hypercube> hyper (new hypercube(Z, X));
    VECEXT::copyHyper(pef, hyper);

    // create extended data vector
    X = dat->getHyper()->getAxis(2);
    Z = dat->getHyper()->getAxis(1);
    X.n += pef->getHyper()->getAxis(2).n -1;
    Z.n += pef->getHyper()->getAxis(1).n -1;
    Z.o = Z.o -Z.d * _lead;
    std::shared_ptr<float2DReg> dataExt (new float2DReg(Z, X));

    // build an equations selector to discard boundaries in the PEF estimation if required
    {
    std::shared_ptr<float2DReg> dat2 (new float2DReg(dat->getHyper()));
    pef->zero();
    dat2->zero();
    VECEXT::add(pef, 1.0);
    VECEXT::add(dat2, 1.0);
    convolution2D conv(pef, dat2->getHyper(), dataExt->getHyper());
    conv.forward(false, dat2, dataExt);
    }

    float N = pef->getHyper()->getAxis(2).n * pef->getHyper()->getAxis(1).n;
    for (int ix=0; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            if ((*dataExt->_mat)[ix][iz] != N) (*dataExt->_mat)[ix][iz] = 0.0;
            else (*dataExt->_mat)[ix][iz] = 1.0;
        }
    }

    selection2D selectEquation(dataExt);
    dataExt->zero();

    // set the pef to zero except for the leading 1
    pef->zero();
    (*pef->_mat)[0][_lead] = 1;

    // define the data convolution matrix
    convolution2D conv(dat, pef->getHyper(), dataExt->getHyper());

    // build the selection operators to separate PEF known and unknown coefficients
    std::shared_ptr<float2DReg> selector (new float2DReg(pef->getHyper()));
    for (int i=0; i<= _lead + _gap; i++)
        (*selector->_mat)[0][i] = 1;

    selection2D selectKnown(selector);
    selection2D selectUnknown(selector);
    selectUnknown.reverse();

    std::shared_ptr<float2DReg> KnP (new float2DReg(pef->getHyper()));
    selectKnown.forward(false, pef, KnP);
    conv.forward(false, KnP, dataExt);
    dataExt->scale(-1);

    // initialize the CGLS solver
    unsigned int niter = VECEXT::getnx(pef)*VECEXT::getnz(pef) - _lead - _gap -1;
    niter = std::min(niter, niter0);
    cg.setNiter(niter);

    std::clog << "Maximum number of iterations for CG solver set to: " <<niter<< "\n";

    chainOper2D chain(&conv, &selectUnknown);

    // Choose whether or not to discard the equations at the bounds and run the solver to estimate the PEF
    if (crossBounds == true){
        cg.run(&chain, pef, dataExt);
    }
    else{
        chainOper2D chain2(&selectEquation, &chain);
        std::shared_ptr<float2DReg> temp = dataExt->clone();
        selectEquation.forward(false, temp, dataExt); 
        cg.run(&chain2, pef, dataExt);
    }
}



void pef2D::estimate(std::shared_ptr<float2DReg> pef, const std::shared_ptr<float2DReg> dat, 
                        const std::shared_ptr<float2DReg> mask, cgls &cg){

    if (_lead > pef->getHyper()->getAxis(1).n-1)
            throw std::logic_error("The lead index exceeds nz");

    if (_lead + _gap > pef->getHyper()->getAxis(1).n-1) 
        throw std::logic_error("Gap too large or filter too short");

    if ((VECEXT::getnx(pef) > VECEXT::getnx(dat)) || (VECEXT::getnz(pef) > VECEXT::getnz(dat)))
        throw std::logic_error("data is smaller than filter");

    dat->getHyper()->checkSame(mask->getHyper());

    // copy the hypercube from data to the PEF without the size and set originx to zero and originz according to the leading 1
    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);
    X.o = 0;
    Z.o = -Z.d * _lead;
    std::shared_ptr<hypercube> hyper (new hypercube(Z, X));
    VECEXT::copyHyper(pef, hyper);

    // create extended data vector
    X = dat->getHyper()->getAxis(2);
    Z = dat->getHyper()->getAxis(1);
    X.n += pef->getHyper()->getAxis(2).n -1;
    Z.n += pef->getHyper()->getAxis(1).n -1;
    Z.o = Z.o -Z.d * _lead;
    std::shared_ptr<float2DReg> dataExt (new float2DReg(Z, X));

    // build an equations selector to discard boundaries and mask in the PEF estimation
    {
    std::shared_ptr<float2DReg> dat2 (new float2DReg(dat->getHyper()));
    std::shared_ptr<float2DReg> dat3 = dat2->clone();
    pef->zero();
    dat2->zero();
    VECEXT::add(pef, 1.0);
    for (int i=0; i<_lead; i++) {(*pef->_mat)[0][i]=0;}
    VECEXT::add(dat2, 1.0);
    selection2D select(mask);
    select.forward(false, dat2, dat3);
    convolution2D conv(pef, dat2->getHyper(), dataExt->getHyper());
    conv.crossBounds(false);
    conv.forward(false, dat2, dataExt);
    }

    float N = pef->getHyper()->getAxis(2).n * pef->getHyper()->getAxis(1).n;
    N = N - _lead;
    for (int ix=0; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            if ((*dataExt->_mat)[ix][iz] != N) (*dataExt->_mat)[ix][iz] = 0.0;
            else (*dataExt->_mat)[ix][iz] = 1.0;
        }
    }

    selection2D selectEquation(dataExt);
    dataExt->zero();

    // set the pef to zero except for the leading 1
    pef->zero();
    (*pef->_mat)[0][_lead] = 1;

    // define the data convolution matrix
    convolution2D conv(dat, pef->getHyper(), dataExt->getHyper());

    // build the selection operators to separate PEF known and unknown coefficients
    std::shared_ptr<float2DReg> selector (new float2DReg(pef->getHyper()));
    for (int i=0; i<= _lead + _gap; i++)
        (*selector->_mat)[0][i] = 1;

    selection2D selectKnown(selector);
    selection2D selectUnknown(selector);
    selectUnknown.reverse();

    {
    std::shared_ptr<float2DReg> KnP (new float2DReg(pef->getHyper()));
    std::shared_ptr<float2DReg> temp = dataExt->clone();
    selectKnown.forward(false, pef, KnP);
    conv.forward(false, KnP, temp);
    selectEquation.forward(false, temp, dataExt);
    dataExt->scale(-1);
    }

    // define the full operator
    chainOper2D chain(&conv, &selectUnknown);
    chainOper2D chain2(&selectEquation, &chain);

    // initialize the CGLS solver
    unsigned int niter = VECEXT::getnx(pef)*VECEXT::getnz(pef) - _lead - _gap -1;
    niter = std::min(niter, cg.getNiter());
    cg.setNiter(niter);

    std::clog << "Maximum number of iterations for CG solver set to: " <<niter<< "\n";

    // solve M.D.U.p = -M.D.K.p     M mask, D data convolution matrix, U unknown pef coef, K known pef coef
    cg.run(&chain2, pef, dataExt);
}