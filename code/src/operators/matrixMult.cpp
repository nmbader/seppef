#include "matrixMult.h"

void matrixMult::forward(bool add, const std::shared_ptr<float1DReg> mod, std::shared_ptr<float1DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nz_mod = mod->getHyper()->getAxis(1).n;
    unsigned int nz_dat = dat->getHyper()->getAxis(1).n;

    if (add == false) dat->scale(0.);

    for (int iz_dat=0; iz_dat<nz_dat; iz_dat++){
        for (int iz_mod=0; iz_mod<nz_mod; iz_mod++){
            (*dat->_mat)[iz_dat] += (*_mat->_mat)[iz_mod][iz_dat] * (*mod->_mat)[iz_mod];
        }
    }
}

void matrixMult::adjoint(bool add, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nz_mod = mod->getHyper()->getAxis(1).n;
    unsigned int nz_dat = dat->getHyper()->getAxis(1).n;

    if (add == false) mod->scale(0.);

    for (int iz_dat=0; iz_dat<nz_dat; iz_dat++){
        for (int iz_mod=0; iz_mod<nz_mod; iz_mod++){
            (*mod->_mat)[iz_mod] += (*_mat->_mat)[iz_mod][iz_dat] * (*dat->_mat)[iz_dat];
        }
    }
}