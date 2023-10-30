#include "identity2D.h"

void identity2D::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nx_mod = mod->getHyper()->getAxis(2).n;
    unsigned int nz_mod = mod->getHyper()->getAxis(1).n;


    if (add == false) dat->scale(0.);

    for (int ix_mod=0; ix_mod<nx_mod; ix_mod++){
        for (int iz_mod=0; iz_mod<nz_mod; iz_mod++){
            (*dat->_mat)[ix_mod][iz_mod] += (*mod->_mat)[ix_mod][iz_mod];
        }
    }
}

void identity2D::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nx_mod = mod->getHyper()->getAxis(2).n;
    unsigned int nz_mod = mod->getHyper()->getAxis(1).n;


    if (add == false) mod->scale(0.);

    for (int ix_mod=0; ix_mod<nx_mod; ix_mod++){
        for (int iz_mod=0; iz_mod<nz_mod; iz_mod++){
            (*mod->_mat)[ix_mod][iz_mod] += (*dat->_mat)[ix_mod][iz_mod];
        }
    }
}