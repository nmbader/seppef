#include <fftw3.h>
#include <complex>
#include <chrono>
#include <thread>
#include "fkTransform.h"
#include "fxTransform.h"
#include "floatHyperExt.h"
#include "complexHyperExt.h"

void fkTransform::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<complex2DReg> dat) const{
    this->checkDomainRange(mod, dat);

    if (add == false) CVECEXT::scale(dat, 0.0);

    unsigned int nf = CVECEXT::getnz(dat);
    unsigned int nz = VECEXT::getnz(mod);
    unsigned int nx = VECEXT::getnx(mod);

    fftwf_plan p;
    float * in = mod->getVals();
    fftwf_complex * out;
    out=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*nf);
      
    p = fftwf_plan_dft_r2c_2d(nx, nz, in, out, FFTW_ESTIMATE);
    fftwf_execute(p);

    std::complex<float>* ptr = reinterpret_cast<std::complex<float>* >(out);

    int ix=0, iz=0;
    for (int i=0; i<nx*nf; i++){
        ix = i/nf;
        iz = i - ix*nf;
        ix += (nx-1)/2;
        ix = ix - nx*(ix/nx);
        //iz += (nf-1)/2;
        //iz = iz - nf*(iz/nf);
        (*dat->_mat)[ix][iz] += ptr[i];
    }
 
    fftwf_destroy_plan(p);
    fftwf_free(out);
    fftwf_cleanup();
}

void fkTransform::inverse(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<complex2DReg> dat) const{
    this->checkDomainRange(mod, dat);

    if (add == false) mod->scale(0.0);

    unsigned int nf = CVECEXT::getnz(dat);
    unsigned int nz = VECEXT::getnz(mod);
    unsigned int nx = VECEXT::getnx(mod);

    fftwf_plan p;
    fftwf_complex * in;
    in=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*nf);
    float * out = (float*)fftwf_malloc(sizeof(float)*nx*nz);
    
    int ix=0, iz=0;
    for (int i=0; i<nx*nf; i++){
        ix = i/nf;
        iz = i - ix*nf;
        ix += (nx-1)/2;
        ix = ix - nx*(ix/nx);
        //iz += (nf-1)/2;
        //iz = iz - nf*(iz/nf);
        in[i][0] = (*dat->_mat)[ix][iz].real();
        in[i][1] = (*dat->_mat)[ix][iz].imag();
    }

    p = fftwf_plan_dft_c2r_2d(nx, nz, in, out, FFTW_ESTIMATE);
    fftwf_execute(p);

    for (int i=0; i<nx*nz; i++){
        ix = i/nz;
        iz = i - ix*nz;
        (*mod->_mat)[ix][iz] = out[i]/(nx*nz);
    }  

    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
    fftwf_cleanup();
}

void fkTransform::cforward(bool add, const std::shared_ptr<complex2DReg> mod, std::shared_ptr<complex2DReg> dat) const{
    this->checkDomainRange(mod, dat);

    if (add == false) CVECEXT::scale(dat, 0.0);

    unsigned int nz = CVECEXT::getnz(mod);
    unsigned int nx = CVECEXT::getnx(mod);

    fftwf_plan p;
    fftwf_complex * in = reinterpret_cast<fftwf_complex* >(mod->getVals());
    fftwf_complex * out;
    out=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*nz);
      
    p = fftwf_plan_dft_2d(nx, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(p);

    std::complex<float>* ptr = reinterpret_cast<std::complex<float>* >(out);

    int ix=0, iz=0;
    for (int i=0; i<nx*nz; i++){
        ix = i/nz;
        iz = i - ix*nz;
        ix += (nx-1)/2;
        ix = ix - nx*(ix/nx);
        iz += (nz-1)/2;
        iz = iz - nz*(iz/nz);
        (*dat->_mat)[ix][iz] += ptr[i];
    }
 
    fftwf_destroy_plan(p);
    fftwf_free(out);
    fftwf_cleanup();
}

void fkTransform::cinverse(bool add, std::shared_ptr<complex2DReg> mod, const std::shared_ptr<complex2DReg> dat) const{
    this->checkDomainRange(mod, dat);

    if (add == false) CVECEXT::scale(mod, 0.0);

    unsigned int nz = CVECEXT::getnz(dat);
    unsigned int nx = CVECEXT::getnx(dat);

    fftwf_plan p;
    fftwf_complex * in;
    fftwf_complex * out;
    in=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*nz);
    out=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nx*nz);

    int ix=0, iz=0;
    for (int i=0; i<nx*nz; i++){
        ix = i/nz;
        iz = i - ix*nz;
        ix += (nx-1)/2;
        ix = ix - nx*(ix/nx);
        iz += (nz-1)/2;
        iz = iz - nz*(iz/nz);
        in[i][0] = (*dat->_mat)[ix][iz].real();
        in[i][1] = (*dat->_mat)[ix][iz].imag();
    }

    p = fftwf_plan_dft_2d(nx, nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(p);

    std::complex<float>* ptr = reinterpret_cast<std::complex<float>* >(out);

    for (int i=0; i<nx*nz; i++){
        ix = i/nz;
        iz = i - ix*nz;
        (*mod->_mat)[ix][iz].real(ptr[i].real()/(nx*nz));
        (*mod->_mat)[ix][iz].imag(ptr[i].imag()/(nx*nz));
    }  

    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
    fftwf_cleanup();
}

void fkTransform::testInverse() {

    std::shared_ptr<float2DReg> mod1 (new float2DReg(_domain));
    std::shared_ptr<float2DReg> mod1b (new float2DReg(_domain));
    axis Z = _domain->getAxis(1);
    axis X = _domain->getAxis(2);
    Z.n += 1;
    X.n += 1;
    std::shared_ptr<float2DReg> mod2 (new float2DReg(Z, X));
    std::shared_ptr<float2DReg> mod2b (new float2DReg(Z, X));

    VECEXT::random(mod1);
    VECEXT::random(mod2);
    
    std::shared_ptr<complex2DReg> dat1 (new complex2DReg(_range));

    this->forward(false,mod1,dat1);
    this->inverse(false,mod1b,dat1);

    this->setDomainRange(mod2);
    std::shared_ptr<complex2DReg> dat2 (new complex2DReg(_range));
    this->forward(false,mod2,dat2);
    this->inverse(false,mod2b,dat2);

    double trueNorm1, trueNorm2, error1, error2;
    trueNorm1 = VECEXT::norm2(mod1);
    trueNorm2 = VECEXT::norm2(mod2);
    mod1b->scale(-1); mod2b->scale(-1);
    mod1b->add(mod1); mod2b->add(mod2);
    error1 = VECEXT::norm2(mod1b);
    error2 = VECEXT::norm2(mod2b);

    std::clog << "relative error for the original parity: " << error1/trueNorm1 << std::endl;
    std::clog << "relative error for the inverse parity: " << error2/trueNorm2 << std::endl;
}

void fkTransform::testCInverse() {

    std::shared_ptr<complex2DReg> mod1 (new complex2DReg(_domain));
    std::shared_ptr<complex2DReg> mod1b (new complex2DReg(_domain));

    CVECEXT::random(mod1);
    
    std::shared_ptr<complex2DReg> dat1 (new complex2DReg(_crange));

    this->cforward(false,mod1,dat1);
    this->cinverse(false,mod1b,dat1);

    double trueNorm1, estimatedNorm1, error1;
    trueNorm1 = CVECEXT::norm2(mod1);
    estimatedNorm1 = CVECEXT::norm2(mod1b);
    CVECEXT::scale(mod1b, -1.0);
    CVECEXT::add(mod1b, mod1);
    error1 = CVECEXT::norm2(mod1b);

    std::clog << "true vector norm2: " << trueNorm1 << std::endl;
    std::clog << "inverted vector norm2: " << estimatedNorm1 << std::endl;
    std::clog << "relative error: " << error1/trueNorm1 << std::endl;
}