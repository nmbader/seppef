#include "identity1D.h"

void identity1D::forward(bool add, const std::shared_ptr<float1DReg> mod, std::shared_ptr<float1DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nz_mod = mod->getHyper()->getAxis(1).n;


    if (add == false) dat->scale(0.);
    for (int iz_mod=0; iz_mod<nz_mod; iz_mod++){
        (*dat->_mat)[iz_mod] += (*mod->_mat)[iz_mod];
    }
}

void identity1D::adjoint(bool add, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nz_mod = mod->getHyper()->getAxis(1).n;


    if (add == false) mod->scale(0.);
    for (int iz_mod=0; iz_mod<nz_mod; iz_mod++){
        (*mod->_mat)[iz_mod] += (*dat->_mat)[iz_mod];
    }
}