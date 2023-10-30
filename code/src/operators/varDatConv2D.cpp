#include "varDatConv2D.h"
#include <omp.h>


void varDatConv2D::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
  
    this->checkDomainRange(mod,dat);

    int nx_f = _fil0->getAxis(2).n;
    int nz_f = _fil0->getAxis(1).n;
    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nz_dat0 = _dat0->getHyper()->getAxis(1).n;

    if (add == false) dat->scale(0.0);

    int nx_mod_left, nz_mod_top;
    float dx, dz;

    float * pdat = dat->getVals();
    float * pmod = mod->getVals();
    float * pdat0 = _dat0->getVals();

    #pragma omp parallel for private(nx_mod_left, nz_mod_top, dx, dz)
    for (int ix_dat= nx_f - 1; ix_dat < nx_dat; ix_dat++){

        nx_mod_left = floor(ix_dat / _xinc);
        dx = (ix_dat - nx_mod_left * _xinc) / (float)_xinc;
        nx_mod_left *= nx_f;

        for (int iz_dat = nz_f -1 -_lead; iz_dat < nz_dat - _lead; iz_dat++){

            nz_mod_top = floor((iz_dat + _lead) / _zinc);
            dz = (iz_dat + _lead - nz_mod_top * _zinc) / (float)_zinc;
            nz_mod_top *= nz_f;

            for (int ix_f = 0; ix_f < nx_f; ix_f++){
                for (int iz_f = 0; iz_f < nz_f; iz_f++){
                    
                    pdat[ix_dat*nz_dat+iz_dat] += pdat0[(ix_dat - ix_f)*nz_dat0+iz_dat - iz_f + _lead]*
                                                    ((1-dx)*(1-dz)*pmod[(nx_mod_left+ix_f)*nz_mod+nz_mod_top+iz_f]+ // upper left filter
                                                    (dz-dx*dz)*pmod[(nx_mod_left+ix_f)*nz_mod+nz_mod_top+nz_f+iz_f]+ // bottom left filter
                                                    (dx-dx*dz)*pmod[(nx_mod_left+nx_f+ix_f)*nz_mod+nz_mod_top+iz_f]+ // upper right filter
                                                    dx*dz*pmod[(nx_mod_left+nx_f+ix_f)*nz_mod+nz_mod_top+nz_f+iz_f]); // bottom right filter

                }
            }
        }
    }
}

void varDatConv2D::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {

    this->checkDomainRange(mod,dat);

    int nx_f = _fil0->getAxis(2).n;
    int nz_f = _fil0->getAxis(1).n;
    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nz_dat0 = _dat0->getHyper()->getAxis(1).n;

    if (add == false) mod->scale(0.0);

    float dx, dz;
    int nx_mod_left, nz_mod_top;

    float * pdat = dat->getVals();
    float * pmod = mod->getVals();
    float * pdat0 = _dat0->getVals();

    for (int ix_dat= nx_f - 1; ix_dat < nx_dat; ix_dat++){

        nx_mod_left = floor(ix_dat / _xinc);
        dx = (ix_dat - nx_mod_left * _xinc) / (float)_xinc;
        nx_mod_left *= nx_f;

        for (int iz_dat = nz_f -1 -_lead; iz_dat < nz_dat - _lead; iz_dat++){

            nz_mod_top = floor((iz_dat + _lead) / _zinc);
            dz = (iz_dat + _lead - nz_mod_top * _zinc) / (float)_zinc;
            nz_mod_top *= nz_f;

            for (int ix_f = 0; ix_f < nx_f; ix_f++){
                for (int iz_f = 0; iz_f < nz_f; iz_f++){

                    pmod[(nx_mod_left+ix_f)*nz_mod+nz_mod_top+iz_f] += pdat0[(ix_dat - ix_f)*nz_dat0+iz_dat - iz_f + _lead]*
                                                                        (1-dx)*(1-dz)*pdat[ix_dat*nz_dat+iz_dat];
                    pmod[(nx_mod_left+ix_f)*nz_mod+nz_mod_top+nz_f+iz_f] += pdat0[(ix_dat - ix_f)*nz_dat0+iz_dat - iz_f + _lead]*
                                                                        (dz-dx*dz)*pdat[ix_dat*nz_dat+iz_dat];
                    pmod[(nx_mod_left+nx_f+ix_f)*nz_mod+nz_mod_top+iz_f] += pdat0[(ix_dat - ix_f)*nz_dat0+iz_dat - iz_f + _lead]*
                                                                        (dx-dx*dz)*pdat[ix_dat*nz_dat+iz_dat];
                    pmod[(nx_mod_left+nx_f+ix_f)*nz_mod+nz_mod_top+nz_f+iz_f] += pdat0[(ix_dat - ix_f)*nz_dat0+iz_dat - iz_f + _lead]*
                                                                        dx*dz*pdat[ix_dat*nz_dat+iz_dat];                                 
                }
            }
        }
    }
}

