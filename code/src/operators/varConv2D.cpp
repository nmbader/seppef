#include "varConv2D.h"

void varConv2D::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nz_f = _allfil->getHyper()->getAxis(1).n;

    if (add == false) dat->scale(0.0);

    int nx_pef_left, nz_pef_top;
    float dx, dz;

    float * pdat = dat->getVals();
    const float * pmod = mod->getCVals();
    const float * pfil = _allfil->getCVals();

    #pragma omp parallel for private(nx_pef_left, nz_pef_top, dx, dz)
    for (int ix_dat= _nx - 1; ix_dat < nx_dat; ix_dat++){

        nx_pef_left = floor(ix_dat / _xinc);
        dx = (ix_dat - nx_pef_left * _xinc) / (float)_xinc;
        nx_pef_left *= _nx;

        for (int iz_dat = _nz -1 -_lead; iz_dat < nz_dat - _lead; iz_dat++){

            nz_pef_top = floor((iz_dat + _lead) / _zinc);
            dz = (iz_dat + _lead - nz_pef_top * _zinc) / (float)_zinc;
            nz_pef_top *= _nz;

            for (int ix_f = 0; ix_f < _nx; ix_f++){
                for (int iz_f = 0; iz_f < _nz; iz_f++){

                    pdat[ix_dat*nz_dat+iz_dat] += pmod[(ix_dat - ix_f)*nz_mod+iz_dat - iz_f + _lead]*
                                                    ((1-dx)*(1-dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper left filter
                                                    (dz-dx*dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+_nz+iz_f]+ // bottom left filter
                                                    (dx-dx*dz)*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper right filter
                                                    dx*dz*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+_nz+iz_f]); // bottom right filter

                }
            }
        }
    }
}

void varConv2D::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nz_f = _allfil->getHyper()->getAxis(1).n;

    if (add == false) mod->scale(0.0);

    int nx_pef_left, nz_pef_top;
    float dx, dz;

    const float * pdat = dat->getCVals();
    float * pmod = mod->getVals();
    const float * pfil = _allfil->getCVals();

    for (int ix_dat= _nx - 1; ix_dat < nx_dat; ix_dat++){

        nx_pef_left = floor(ix_dat / _xinc);
        dx = (ix_dat - nx_pef_left * _xinc) / (float)_xinc;
        nx_pef_left *= _nx;

        for (int iz_dat = _nz -1 -_lead; iz_dat < nz_dat - _lead; iz_dat++){

            nz_pef_top = floor((iz_dat + _lead) / _zinc);
            dz = (iz_dat + _lead - nz_pef_top * _zinc) / (float)_zinc;
            nz_pef_top *= _nz;

            for (int ix_f = 0; ix_f < _nx; ix_f++){
                for (int iz_f = 0; iz_f < _nz; iz_f++){

                    pmod[(ix_dat - ix_f)*nz_mod+iz_dat - iz_f + _lead] += pdat[ix_dat*nz_dat+iz_dat]*
                                                    ((1-dx)*(1-dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper left filter
                                                    (dz-dx*dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+_nz+iz_f]+ // bottom left filter
                                                    (dx-dx*dz)*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper right filter
                                                    dx*dz*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+_nz+iz_f]); // bottom right filter

                }
            }
        }
    }
}