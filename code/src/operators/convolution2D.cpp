#include "convolution2D.h"
#include <omp.h>


void convolution2D::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    int nx_f = _filter->getHyper()->getAxis(2).n;
    int nz_f = _filter->getHyper()->getAxis(1).n;
    int nx_mod = mod->getHyper()->getAxis(2).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;

    if (add == false) dat->scale(0.0);

    int max_nx_f, max_nx_mod, max_nz_f, max_nz_mod;
    max_nx_f = std::min(nx_f, nx_dat);
    max_nz_f = std::min(nz_f, nz_dat);

    float * pdat = dat->getVals();
    const float * pmod = mod->getCVals();
    const float * pfil = _filter->getCVals();

    if (_crossBounds == true){
    // Rolling the filter over the interior of the model and the bounds
        for (int ix_f=0; ix_f<max_nx_f; ix_f++){
            max_nx_mod = std::min(nx_dat - ix_f, nx_mod);
            for (int ix_mod=0; ix_mod<max_nx_mod; ix_mod++){
                for (int iz_f=0; iz_f<max_nz_f; iz_f++){
                    max_nz_mod = std::min(nz_dat - iz_f, nz_mod);
                    for (int iz_mod=0; iz_mod<max_nz_mod; iz_mod++){
                        pdat[(ix_f+ix_mod)*nz_dat+iz_f+iz_mod] += pfil[ix_f*nz_f+iz_f] * pmod[ix_mod*nz_mod+iz_mod];
                    }
                }
            }
        }
    }

    else {
    // Rolling the filter over the interior of the model
        int min_nx_f, min_nz_f;
        max_nx_mod = std::min(nx_mod, nx_dat);
        max_nz_mod = std::min(nz_mod, nz_dat);
        for (int ix_mod=0; ix_mod<max_nx_mod; ix_mod++){
            min_nx_f = std::max(0, nx_f - 1 -ix_mod);
            max_nx_f = std::min(nx_f, max_nx_mod - ix_mod);
            for (int ix_f=min_nx_f; ix_f<max_nx_f; ix_f++){
                for (int iz_mod=0; iz_mod<max_nz_mod; iz_mod++){
                    min_nz_f = std::max(0, nz_f - 1 -iz_mod);
                    max_nz_f = std::min(nz_f, max_nz_mod - iz_mod);
                    for (int iz_f=min_nz_f; iz_f<max_nz_f; iz_f++){
                        pdat[(ix_f+ix_mod)*nz_dat+iz_f+iz_mod] += pfil[ix_f*nz_f+iz_f] * pmod[ix_mod*nz_mod+iz_mod];
                    }
                }
            }
        }
    }
}

void convolution2D::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    int nx_f = _filter->getHyper()->getAxis(2).n;
    int nz_f = _filter->getHyper()->getAxis(1).n;
    int nx_mod = mod->getHyper()->getAxis(2).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;

    if (add == false) mod->scale(0.0);

    int max_nx_f, max_nx_mod, max_nz_f, max_nz_mod;
    max_nx_f = std::min(nx_f, nx_dat);
    max_nz_f = std::min(nz_f, nz_dat);

    // get the total number of threads
    int nth = omp_get_num_threads();
    int nth1 = nth/2 + 1;
    int nth2 = nth - nth1;

    const float * pdat = dat->getCVals();
    float * pmod = mod->getVals();
    const float * pfil = _filter->getCVals();

    if (_crossBounds == true){
    // Rolling the filter over the interior of the model and the bounds
        for (int ix_f=0; ix_f<max_nx_f; ix_f++){
            max_nx_mod = std::min(nx_dat - ix_f, nx_mod);
            #pragma omp parallel for private(max_nz_mod)
            for (int ix_mod=0; ix_mod<max_nx_mod; ix_mod++){
                for (int iz_f=0; iz_f<max_nz_f; iz_f++){
                    max_nz_mod = std::min(nz_dat - iz_f, nz_mod);
                    for (int iz_mod=0; iz_mod<max_nz_mod; iz_mod++){
                        pmod[ix_mod*nz_mod+iz_mod] += pfil[ix_f*nz_f+iz_f] * pdat[(ix_f+ix_mod)*nz_dat+iz_f+iz_mod];
                    }
                }
            }
        }
    }

    else {
    // Rolling the filter over the interior of the model
        int min_nx_f, min_nz_f;
        max_nx_mod = std::min(nx_mod, nx_dat);
        max_nz_mod = std::min(nz_mod, nz_dat);
        #pragma omp parallel for private(min_nx_f, max_nx_f, min_nz_f, max_nz_f)
        for (int ix_mod=0; ix_mod<max_nx_mod; ix_mod++){
            min_nx_f = std::max(0, nx_f - 1 -ix_mod);
            max_nx_f = std::min(nx_f, max_nx_mod - ix_mod);
            for (int ix_f=min_nx_f; ix_f<max_nx_f; ix_f++){
                for (int iz_mod=0; iz_mod<max_nz_mod; iz_mod++){
                    min_nz_f = std::max(0, nz_f - 1 -iz_mod);
                    max_nz_f = std::min(nz_f, max_nz_mod - iz_mod);
                    for (int iz_f=min_nz_f; iz_f<max_nz_f; iz_f++){
                        pmod[ix_mod*nz_mod+iz_mod] += pfil[ix_f*nz_f+iz_f] * pdat[(ix_f+ix_mod)*nz_dat+iz_f+iz_mod];
                    }
                }
            }
        }
    }
}