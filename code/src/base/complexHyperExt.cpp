
#include "complexHyperExt.h"
#include "floatHyperExt.h"
#include <stdexcept>
#include <cstdlib>


using namespace SEP;

unsigned int CVECEXT::getnx(const std::shared_ptr<complex2DReg> cvec) {
    return cvec->getHyper()->getAxis(2).n;}

unsigned int CVECEXT::getnz(const std::shared_ptr<complex2DReg> cvec) {
    return cvec->getHyper()->getAxis(1).n;}

std::shared_ptr<complex2DReg> CVECEXT::makeComplex (const std::shared_ptr<float2DReg> vec){
    std::shared_ptr<complex2DReg> cvec (new complex2DReg(vec->getHyper()));
    unsigned int nx = VECEXT::getnx(vec);
    unsigned int nz = VECEXT::getnz(vec);
    std::shared_ptr<float2D> vals = vec->_mat;
    std::shared_ptr<complex2D> cvals = cvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz; iz++){
            (*cvals)[ix][iz].real((*vals)[ix][iz]);
            (*cvals)[ix][iz].imag(0);
        }
    }
    return cvec;
}

void CVECEXT::add(std::shared_ptr<complex2DReg> cvec, const float z){
    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals)[ix][iz].real(z+(*cvals)[ix][iz].real());
        }
    }    
}

void CVECEXT::scale(std::shared_ptr<complex2DReg> cvec, const float z){
    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;   
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals)[ix][iz].real(z*(*cvals)[ix][iz].real());
            (*cvals)[ix][iz].imag(z*(*cvals)[ix][iz].imag());
        }
    }    
}

void CVECEXT::add(std::shared_ptr<complex2DReg> cvec, const std::complex<float> z){
    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;    
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals)[ix][iz] += z;
        }
    }    
}

void CVECEXT::scale(std::shared_ptr<complex2DReg> cvec, const std::complex<float> z){
    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;    
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals)[ix][iz] *= z;
        }
    }    
}

void CVECEXT::add(std::shared_ptr<complex2DReg> cvec1, const std::shared_ptr<complex2DReg> cvec2){
    cvec1->checkSame(cvec2);
    unsigned int nx = CVECEXT::getnx(cvec1);
    unsigned int nz = CVECEXT::getnz(cvec1);
    std::shared_ptr<complex2D> cvals1 = cvec1->_mat;
    std::shared_ptr<complex2D> cvals2 = cvec2->_mat;
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals1)[ix][iz] = (*cvals1)[ix][iz] + (*cvals2)[ix][iz];
        }
    }  
}

void CVECEXT::scaleAdd(std::shared_ptr<complex2DReg> cvec1, const std::shared_ptr<complex2DReg> cvec2, const std::complex<float> alpha, const std::complex<float> beta){
    cvec1->checkSame(cvec2);
    unsigned int nx = CVECEXT::getnx(cvec1);
    unsigned int nz = CVECEXT::getnz(cvec1);
    std::shared_ptr<complex2D> cvals1 = cvec1->_mat;
    std::shared_ptr<complex2D> cvals2 = cvec2->_mat;
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals1)[ix][iz] = alpha*(*cvals1)[ix][iz] + beta*(*cvals2)[ix][iz];
        }
    }    
}

void CVECEXT::random(std::shared_ptr<complex2DReg> cvec, float min, float max){
    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals)[ix][iz].real(min + (max-min)*(rand() % 10000)/10000.0);
            (*cvals)[ix][iz].imag(min + (max-min)*(rand() % 10000)/10000.0);
        }
    }    
}

std::complex<double> CVECEXT::dot(const std::shared_ptr<complex2DReg> cvec1, const std::shared_ptr<complex2DReg> cvec2){
    cvec1->checkSame(cvec2);
    unsigned int nx = CVECEXT::getnx(cvec1);
    unsigned int nz = CVECEXT::getnz(cvec1);
    std::shared_ptr<complex2D> cvals1 = cvec1->_mat;
    std::shared_ptr<complex2D> cvals2 = cvec2->_mat;
    std::complex<double> sum = 0;  
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            sum += std::conj((*cvals1)[ix][iz]) * (*cvals2)[ix][iz];
        }
    }
    return sum;    
}

double CVECEXT::norm2(const std::shared_ptr<complex2DReg> cvec){
    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    double sum = 0;  
    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            sum += std::norm((*cvals)[ix][iz]);
        }
    }
    return sum; 
}

void CVECEXT::conjugate(std::shared_ptr<complex2DReg> cvec){
    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*cvals)[ix][iz] = std::conj((*cvals)[ix][iz]);
        }
    }
}

std::shared_ptr<complex2DReg> CVECEXT::adjoint(const std::shared_ptr<complex2DReg> cvec){
    axis X = cvec->getHyper()->getAxis(2);
    axis Z = cvec->getHyper()->getAxis(1);
    std::shared_ptr<complex2DReg> newcvec (new complex2DReg(X, Z));

    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    std::shared_ptr<complex2D> newcvals = newcvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*newcvals)[iz][ix] = std::conj((*cvals)[ix][iz]);
        }
    }
    return newcvec;
}

std::shared_ptr<float2DReg> CVECEXT::real(const std::shared_ptr<complex2DReg> cvec){
    axis X = cvec->getHyper()->getAxis(2);
    axis Z = cvec->getHyper()->getAxis(1);
    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));

    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*newvals)[ix][iz]= (*cvals)[ix][iz].real();
        }
    }
    return newvec;
}

std::shared_ptr<float2DReg> CVECEXT::imag(const std::shared_ptr<complex2DReg> cvec){
    axis X = cvec->getHyper()->getAxis(2);
    axis Z = cvec->getHyper()->getAxis(1);
    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));

    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*newvals)[ix][iz]= (*cvals)[ix][iz].imag();
        }
    }
    return newvec;
}

std::shared_ptr<float2DReg> CVECEXT::module(const std::shared_ptr<complex2DReg> cvec){
    axis X = cvec->getHyper()->getAxis(2);
    axis Z = cvec->getHyper()->getAxis(1);
    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));

    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*newvals)[ix][iz]= std::abs((*cvals)[ix][iz]);
        }
    }
    return newvec;
}

std::shared_ptr<float2DReg> CVECEXT::module2(const std::shared_ptr<complex2DReg> cvec){
    axis X = cvec->getHyper()->getAxis(2);
    axis Z = cvec->getHyper()->getAxis(1);
    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));

    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*newvals)[ix][iz]= std::norm((*cvals)[ix][iz]);
        }
    }
    return newvec;
}

std::shared_ptr<float2DReg> CVECEXT::phase(const std::shared_ptr<complex2DReg> cvec){
    axis X = cvec->getHyper()->getAxis(2);
    axis Z = cvec->getHyper()->getAxis(1);
    std::shared_ptr<float2DReg> newvec (new float2DReg(Z, X));

    unsigned int nx = CVECEXT::getnx(cvec);
    unsigned int nz = CVECEXT::getnz(cvec);
    std::shared_ptr<complex2D> cvals = cvec->_mat;
    std::shared_ptr<float2D> newvals = newvec->_mat;

    for (int ix=0; ix<nx; ix++){
        for (int iz=0; iz<nz ; iz++){
            (*newvals)[ix][iz]= std::arg((*cvals)[ix][iz]);
        }
    }
    return newvec;
}

std::shared_ptr<complex2DReg> CVECEXT::resize(const std::shared_ptr<complex2DReg> vec, int nx, int nz){
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    X.n = nx;
    Z.n = nz;
    std::shared_ptr<complex2DReg> newvec (new complex2DReg(Z, X));

    std::shared_ptr<complex2D> vals = vec->_mat;
    std::shared_ptr<complex2D> newvals = newvec->_mat;

    for (int ix=0; ix < std::min(nx,vec->getHyper()->getAxis(2).n); ix++){
        for (int iz=0; iz < std::min(nz,vec->getHyper()->getAxis(1).n); iz++){
            (*newvals)[ix][iz] = (*vals)[ix][iz];
        }
    }
    return newvec;
}

std::shared_ptr<complex2DReg> CVECEXT::pad(const std::shared_ptr<complex2DReg> vec, int nx_start, int nx_end,
                                unsigned int nz_start, unsigned int nz_end){
    axis X = vec->getHyper()->getAxis(2);
    axis Z = vec->getHyper()->getAxis(1);
    unsigned int nx = X.n;
    unsigned int nz = Z.n;
    X.n += nx_start + nx_end;
    Z.n += nz_start + nz_end;
    X.o = X.o - nx_start*X.d;
    Z.o = Z.o - nz_start*Z.d;

    std::shared_ptr<complex2DReg> newvec (new complex2DReg(Z, X));
    newvec->zero();

    std::shared_ptr<complex2D> vals = vec->_mat;
    std::shared_ptr<complex2D> newvals = newvec->_mat;

    for (int ix=nx_start; ix<nx_start+nx; ix++){
        for (int iz=nz_start; iz<nz_start+nz; iz++){
            (*newvals)[ix][iz] = (*vals)[ix-nx_start][iz-nz_start];
        }
    }
    return newvec;
}

void CVECEXT::copyHyper(std::shared_ptr<complex2DReg> vec, const std::shared_ptr<hypercube> hyper){
    std::vector<axis> axes = hyper->getAxes();
    axes[1].n = vec->getHyper()->getAxis(2).n;
    axes[0].n = vec->getHyper()->getAxis(1).n;

    std::shared_ptr<hypercube> h (new hypercube(axes[0], axes[1]));
    vec->setHyper(h);
}