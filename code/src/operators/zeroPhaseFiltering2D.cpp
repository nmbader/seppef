#include "zeroPhaseFiltering2D.h"


void zeroPhaseFiltering2D::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nz_f = _filter->getHyper()->getAxis(1).n/2;
    unsigned int nx = dat->getHyper()->getAxis(2).n;
    unsigned int nz = dat->getHyper()->getAxis(1).n;

    if (add == false) dat->scale(0.);

    for (int ix=0; ix < nx; ix++){

    // apply filtering in the roll-on part
        for (int iz=nz_f; iz < 2*nz_f+1; iz++){
            for (int iz_f=0; iz_f<=iz; iz_f++){
                (*dat->_mat)[ix][iz-nz_f] += (*mod->_mat)[ix][iz-iz_f]*(*_filter->_mat)[iz_f];
            }
        }

    // apply filtering to the rest of the data
        for (int iz=2*nz_f+1; iz < nz; iz++){
            for (int iz_f=0; iz_f<2*nz_f+1; iz_f++){
                (*dat->_mat)[ix][iz-nz_f] += (*mod->_mat)[ix][iz-iz_f]*(*_filter->_mat)[iz_f];
            }
        }

    // apply filtering in the roll-off part
        for (int iz=nz; iz < nz+nz_f; iz++){
            for (int iz_f=iz-nz+1; iz_f<2*nz_f+1; iz_f++){
                (*dat->_mat)[ix][iz-nz_f] += (*mod->_mat)[ix][iz-iz_f]*(*_filter->_mat)[iz_f];
            }
        }

    }

}

void zeroPhaseFiltering2D::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    unsigned int nz_f = _filter->getHyper()->getAxis(1).n/2;
    unsigned int nx = dat->getHyper()->getAxis(2).n;
    unsigned int nz = dat->getHyper()->getAxis(1).n;

    if (add == false) mod->scale(0.);

    for (int ix=0; ix < nx; ix++){

    // apply adjoint filtering in the roll-on part
        for (int iz=nz_f; iz < 2*nz_f+1; iz++){
            for (int iz_f=0; iz_f<=iz; iz_f++){
                (*mod->_mat)[ix][iz-iz_f] += (*dat->_mat)[ix][iz-nz_f]*(*_filter->_mat)[iz_f];
            }
        }

    // apply adjoint filtering to the rest of the data
        for (int iz=2*nz_f+1; iz < nz; iz++){
            for (int iz_f=0; iz_f<2*nz_f+1; iz_f++){
                (*mod->_mat)[ix][iz-iz_f] += (*dat->_mat)[ix][iz-nz_f]*(*_filter->_mat)[iz_f];
            }
        }

    // apply adjoint filtering in the roll-off part
        for (int iz=nz; iz < nz+nz_f; iz++){
            for (int iz_f=iz-nz+1; iz_f<2*nz_f+1; iz_f++){
                (*mod->_mat)[ix][iz-iz_f] += (*dat->_mat)[ix][iz-nz_f]*(*_filter->_mat)[iz_f];
            }
        }

    }

}
