#include "expVarConv2D.h"

void expVarConv2D::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nz_f = _allfil->getHyper()->getAxis(1).n;

    if (add == false) dat->scale(0.0);

    int nx_pef_left, nz_pef_top;
    int nx_f_max, nz_f_max, nz_f_min;
    float dx, dz;

    float * pdat = dat->getVals();
    const float * pmod = mod->getCVals();
    const float * pfil = _allfil->getCVals();

    #pragma omp parallel for private(nx_pef_left, nz_pef_top, dx, dz, nx_f_max, nz_f_max, nz_f_min)
    for (int ix_dat= _nx - 1; ix_dat < nx_dat; ix_dat++){

        nx_pef_left = floor(ix_dat / _xinc);
        dx = (ix_dat - nx_pef_left * _xinc) / (float)_xinc;
        nx_pef_left *= _nx;

        for (int iz_dat = _nz -1 -_lead; iz_dat < nz_dat - _lead; iz_dat++){

            nz_pef_top = floor((iz_dat + _lead) / _zinc);
            dz = (iz_dat + _lead - nz_pef_top * _zinc) / (float)_zinc;
            nz_pef_top *= _nz;

            nx_f_max=std::min(_nx, ix_dat/2 + 1);
            nz_f_max=std::min(_nz, 1 + (iz_dat+2*(int)_lead)/2);
            nz_f_min=std::max(0, (iz_dat+2*(int)_lead-nz_dat+2)/2);

            for (int ix_f = 0; ix_f < nx_f_max; ix_f++){
                for (int iz_f = nz_f_min; iz_f < nz_f_max; iz_f++){

                    pdat[ix_dat*nz_dat+iz_dat] += pmod[(ix_dat - 2*ix_f)*nz_mod+iz_dat - 2*iz_f + 2*_lead]*
                                                    ((1-dx)*(1-dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper left filter
                                                    (dz-dx*dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+_nz+iz_f]+ // bottom left filter
                                                    (dx-dx*dz)*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper right filter
                                                    dx*dz*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+_nz+iz_f]); // bottom right filter

                }
            }
        }
    }
}

void expVarConv2D::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    int nx_dat = dat->getHyper()->getAxis(2).n;
    int nz_dat = dat->getHyper()->getAxis(1).n;
    int nz_mod = mod->getHyper()->getAxis(1).n;
    int nz_f = _allfil->getHyper()->getAxis(1).n;

    if (add == false) mod->scale(0.0);

    int nx_pef_left, nz_pef_top;
    int nx_f_max, nz_f_max, nz_f_min;
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

            nx_f_max=std::min(_nx, ix_dat/2 + 1);
            nz_f_max=std::min(_nz, 1 + (iz_dat+2*(int)_lead)/2);
            nz_f_min=std::max(0, (iz_dat+2*(int)_lead-nz_dat+2)/2);

            for (int ix_f = 0; ix_f < nx_f_max; ix_f++){
                for (int iz_f = nz_f_min; iz_f < nz_f_max; iz_f++){

                    pmod[(ix_dat - 2*ix_f)*nz_mod+iz_dat - 2*iz_f + 2*_lead] += pdat[ix_dat*nz_dat+iz_dat]*
                                                    ((1-dx)*(1-dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper left filter
                                                    (dz-dx*dz)*pfil[(nx_pef_left+ix_f)*nz_f+nz_pef_top+_nz+iz_f]+ // bottom left filter
                                                    (dx-dx*dz)*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+iz_f]+ // upper right filter
                                                    dx*dz*pfil[(nx_pef_left+_nx+ix_f)*nz_f+nz_pef_top+_nz+iz_f]); // bottom right filter

                }
            }
        }
    }
}