#include "nsPef2D.h"
#include "varDatConv2D.h"
#include "selection2D.h"
#include "chainOper2D.h"

std::shared_ptr<float2DReg> nsPef2D_v0::computeDatNorm2(){

    std::shared_ptr<float2DReg> Norm2 = _dat0->clone();
    Norm2->zero();

    // define arrays size
    unsigned int nx_p = _pef->getHyper()->getAxis(2).n;
    unsigned int nz_p = _pef->getHyper()->getAxis(1).n;
    unsigned int nx_d = _dat0->getHyper()->getAxis(2).n;
    unsigned int nz_d = _dat0->getHyper()->getAxis(1).n;

    double norm2=0;

    for (int ix = nx_p -1 ; ix<nx_d; ix++){
        for (int iz=nz_p-_lead-1; iz<nz_d-_lead; iz++){
            norm2 = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    norm2 += (*_dat0->_mat)[ix-ix_p][iz - iz_p + _lead] * (*_dat0->_mat)[ix-ix_p][iz - iz_p +_lead];
                }
            }
            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                norm2 -= (*_dat0->_mat)[ix][iz - iz_p + _lead] * (*_dat0->_mat)[ix][iz - iz_p + _lead];
            }
            (*Norm2->_mat)[ix][iz] = norm2;
        }
    }
    return Norm2;
}

void nsPef2D_v0::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) dat->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double res_pef, step_length;    


    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (*pef->_mat)[ix_p][iz_p] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
            }

            // update the pef
            step_length = _epsilon;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (*pef->_mat)[ix_p][iz_p] - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dat->_mat)[ix_r][iz_r] += (*pef->_mat)[ix_p][iz_p] * (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

        }
    }
    _peff = pef->clone();

}

void nsPef2D_v0::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) mod->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double res_pef, step_length;    

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (*pef->_mat)[ix_p][iz_p] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
            }

            // update the pef
            step_length = _epsilon;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (*pef->_mat)[ix_p][iz_p] - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] += (*pef->_mat)[ix_p][iz_p] * (*dat->_mat)[ix_r][iz_r];
                }
            }

        }
    }
    _peff = pef->clone();
}

void nsPef2D_v1::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) dat->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double res_pef, data_norm, step_length;    

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (*pef->_mat)[ix_p][iz_p] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = res_pef*res_pef/(res_pef*res_pef*data_norm + _epsilon);
            pef->scaleAdd(dpef, 1, -step_length);
            /* for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (*pef->_mat)[ix_p][iz_p] - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            } */

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dat->_mat)[ix_r][iz_r] += (*pef->_mat)[ix_p][iz_p] * (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

        }
    }
    _peff = pef->clone();

}

void nsPef2D_v1::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) mod->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double res_pef, data_norm, step_length;    


    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (*pef->_mat)[ix_p][iz_p] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = res_pef*res_pef/(res_pef*res_pef*data_norm + _epsilon);
            pef->scaleAdd(dpef, 1, -step_length);
            /* for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (*pef->_mat)[ix_p][iz_p] - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            } */

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] += (*pef->_mat)[ix_p][iz_p] * (*dat->_mat)[ix_r][iz_r];
                }
            }

        }
    }
    _peff = pef->clone();
}


void nsPef2D_v2::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) dat->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> pefx;
    std::shared_ptr<float2DReg> pefz;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();

    double res_pef, data_norm, step_length;
    int array_size;

    // store an array of pef of the same size as nz_d
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // erase the first PEF
            if (ix_r > nx_p - 1)
                pef_array.erase(pef_array.begin());

            array_size = pef_array.size();
            //pefz = pef_array[(array_size-1)*std::min(1, iz_r - nz_p + (int)_lead + 1)];
            pefz = pef_array[array_size-1];
            pefx = pef_array[(1-std::min(1, ix_r - nx_p + 1))*(array_size-1)];

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p])
                    * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = res_pef*res_pef/(res_pef*res_pef*data_norm + _epsilon);
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p]) 
                    - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // store the new pef
            newpef = pef->clone();
            pef_array.push_back(newpef);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dat->_mat)[ix_r][iz_r] += (*pef->_mat)[ix_p][iz_p] * (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

        }

    }
    _peff = pef->clone();

}

void nsPef2D_v2::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) mod->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> pefx;
    std::shared_ptr<float2DReg> pefz;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();

    double res_pef, data_norm, step_length;
    int array_size;

    // store an array of pef of the same size as nz_d
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // erase the first PEF
            if (ix_r > nx_p - 1)
                pef_array.erase(pef_array.begin());

            array_size = pef_array.size();
            //pefz = pef_array[(array_size-1)*std::min(1, iz_r - nz_p + (int)_lead + 1)];
            pefz = pef_array[array_size-1];
            pefx = pef_array[(1-std::min(1, ix_r - nx_p + 1))*(array_size-1)];

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p])
                    * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = res_pef*res_pef/(res_pef*res_pef*data_norm + _epsilon);
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p]) 
                    - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // store the new pef
            newpef = pef->clone();
            pef_array.push_back(newpef);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] += (*pef->_mat)[ix_p][iz_p] * (*dat->_mat)[ix_r][iz_r];
                }
            }

        }
    }
    _peff = pef->clone();
}



void nsPef2D_v3::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) dat->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> pefx;
    std::shared_ptr<float2DReg> pefz;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();

    double res_pef, data_norm, step_length;
    int array_size;

    // store an array of pef of the same size as nz_d
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // erase the first PEF
            if (ix_r > nx_p - 1)
                pef_array.erase(pef_array.begin());

            array_size = pef_array.size();
            //pefz = pef_array[(array_size-1)*std::min(1, iz_r - nz_p + (int)_lead + 1)];
            pefz = pef_array[array_size-1];
            pefx = pef_array[(1-std::min(1, ix_r - nx_p + 1))*(array_size-1)];

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p])
                    * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = 1/(_epsilon*data_norm + (1-_epsilon)*_d2max);
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p]) 
                    - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // store the new pef
            newpef = pef->clone();
            pef_array.push_back(newpef);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dat->_mat)[ix_r][iz_r] += (*pef->_mat)[ix_p][iz_p] * (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

        }

    }
    _peff = pef->clone();

}

void nsPef2D_v3::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) mod->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> pefx;
    std::shared_ptr<float2DReg> pefz;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();

    double res_pef, data_norm, step_length;
    int array_size;

    // store an array of pef of the same size as nz_d
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // erase the first PEF
            if (ix_r > nx_p - 1)
                pef_array.erase(pef_array.begin());

            array_size = pef_array.size();
            //pefz = pef_array[(array_size-1)*std::min(1, iz_r - nz_p + (int)_lead + 1)];
            pefz = pef_array[array_size-1];
            pefx = pef_array[(1-std::min(1, ix_r - nx_p + 1))*(array_size-1)];

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p])
                    * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = 1/(_epsilon*data_norm + (1-_epsilon)*_d2max);
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p]) 
                    - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // store the new pef
            newpef = pef->clone();
            pef_array.push_back(newpef);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] += (*pef->_mat)[ix_p][iz_p] * (*dat->_mat)[ix_r][iz_r];
                }
            }

        }
    }
    _peff = pef->clone();
}



void nsPef2D_v4::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) dat->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double res_pef, data_norm, step_length;    

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (*pef->_mat)[ix_p][iz_p] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = 1/(data_norm + _epsilon);
            pef->scaleAdd(dpef, 1, -step_length);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dat->_mat)[ix_r][iz_r] += (*pef->_mat)[ix_p][iz_p] * (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

        }
    }
    _peff = pef->clone();

}

void nsPef2D_v4::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) mod->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double res_pef, data_norm, step_length;    


    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (*pef->_mat)[ix_p][iz_p] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = 1/(data_norm + _epsilon);
            pef->scaleAdd(dpef, 1, -step_length);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] += (*pef->_mat)[ix_p][iz_p] * (*dat->_mat)[ix_r][iz_r];
                }
            }

        }
    }
    _peff = pef->clone();
}

void nsPef2D_v5::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) dat->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> pefx;
    std::shared_ptr<float2DReg> pefz;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();

    double res_pef, data_norm, step_length;
    int array_size;

    // store an array of pef of the same size as nz_d
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // erase the first PEF
            if (ix_r > nx_p - 1)
                pef_array.erase(pef_array.begin());

            array_size = pef_array.size();
            //pefz = pef_array[(array_size-1)*std::min(1, iz_r - nz_p + (int)_lead + 1)];
            pefz = pef_array[array_size-1];
            pefx = pef_array[(1-std::min(1, ix_r - nx_p + 1))*(array_size-1)];

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p])
                    * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = 1/(data_norm + _epsilon/(res_pef*res_pef));
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p]) 
                    - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // store the new pef
            newpef = pef->clone();
            pef_array.push_back(newpef);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dat->_mat)[ix_r][iz_r] += (*pef->_mat)[ix_p][iz_p] * (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

        }

    }
    _peff = pef->clone();

}

void nsPef2D_v5::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) mod->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define intermediate vectors
    std::shared_ptr<float2DReg> pef = _pef->clone();
    std::shared_ptr<float2DReg> pefx;
    std::shared_ptr<float2DReg> pefz;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();

    double res_pef, data_norm, step_length;
    int array_size;

    // store an array of pef of the same size as nz_d
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//

    // loop over the columns
    for (int ix_r=nx_p - 1 ; ix_r<nx_d; ix_r++){

        // loop over the rows
        for (int iz_r=nz_p - _lead - 1; iz_r<nz_d - _lead; iz_r++){

            // erase the first PEF
            if (ix_r > nx_p - 1)
                pef_array.erase(pef_array.begin());

            array_size = pef_array.size();
            //pefz = pef_array[(array_size-1)*std::min(1, iz_r - nz_p + (int)_lead + 1)];
            pefz = pef_array[array_size-1];
            pefx = pef_array[(1-std::min(1, ix_r - nx_p + 1))*(array_size-1)];

            // compute residuals
            res_pef = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    res_pef += (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p])
                    * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            // compute filter perturbation
            dpef->zero();
            data_norm = 0;
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*dpef->_mat)[ix_p][iz_p] += res_pef * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                    data_norm += (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r - ix_p][iz_r - iz_p + _lead];
                }
            }

            for (int iz_p=0; iz_p<_lead+1; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
                data_norm -= (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead] * (*_dat0->_mat)[ix_r][iz_r - iz_p + _lead];
            }

            // update the pef
            step_length = 1/(data_norm + _epsilon/(res_pef*res_pef));
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*pef->_mat)[ix_p][iz_p] = (_weightx*(*pefx->_mat)[ix_p][iz_p]+_weightz*(*pefz->_mat)[ix_p][iz_p]) 
                    - step_length*(*dpef->_mat)[ix_p][iz_p];
                }
            }

            // store the new pef
            newpef = pef->clone();
            pef_array.push_back(newpef);

            // apply the new PEF to the output
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int iz_p=0; iz_p<nz_p; iz_p++){
                    (*mod->_mat)[ix_r - ix_p][iz_r - iz_p + _lead] += (*pef->_mat)[ix_p][iz_p] * (*dat->_mat)[ix_r][iz_r];
                }
            }

        }
    }
    _peff = pef->clone();
}

void nsPef2D_v9::forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) dat->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define a new structure for the local blocks
    struct block {
        int xstart, xend, zstart, zend, xsize, zsize;
    };

    // define the number of local blocks
    int nwinx = (nx_d - _winx) / (_winx - nx_p + 1);
    nwinx++;
    int nwinz = (nz_d - _winz) / (_winz - nz_p + 1);
    nwinz++;

//std::clog << "Number of blocks in x = " << nwinx << "\n";
//std::clog << "Number of blocks in z = " << nwinz << "\n";

    // build the list of local blocks
    std::vector<block> blocks (nwinx*nwinz);
    for (int ix=0; ix<nwinx; ix++){
        for (int iz=0; iz<nwinz; iz++){
            blocks[ix*nwinz+iz].xstart = ix*(_winx - nx_p + 1);
            blocks[ix*nwinz+iz].zstart = iz*(_winz - nz_p + 1);
            blocks[ix*nwinz+iz].xend = blocks[ix*nwinz+iz].xstart + _winx; 
            blocks[ix*nwinz+iz].zend = blocks[ix*nwinz+iz].zstart + _winz;

            if (ix == nwinx - 1){
                blocks[ix*nwinz+iz].xend = nx_d;
                }
            if (iz == nwinz - 1){
                blocks[ix*nwinz+iz].zend = nz_d;
            }

            blocks[ix*nwinz+iz].xsize = blocks[ix*nwinz+iz].xend - blocks[ix*nwinz+iz].xstart;
            blocks[ix*nwinz+iz].zsize = blocks[ix*nwinz+iz].zend - blocks[ix*nwinz+iz].zstart;
        }
    }

//std::clog << "Maximum block size in x = " << blocks[nwinx*nwinz-1].xsize << "\n";
//std::clog << "Maximum block size in z = " << blocks[nwinx*nwinz-1].zsize << "\n";

    // define residual vectors for the largest block (last one)
    std::vector<std::vector <float> > res_temp (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));
    std::vector<std::vector <float> > res_temp2 (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));
    std::vector<std::vector <float> > res_temp3 (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));
    std::vector<std::vector <float> > res_temp4 (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));

    // define intermediate vectors and variables
    std::shared_ptr<float2DReg> pef;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double step_length, norm2_dg, norm2_g, norm2_res, norm2_res2, norm2_res3, norm2_res4, norm2_dat, w_up, w_left, w_up_left, w_dn_left, w_sum, dx, dz;
    int xstart, xend, zstart, zend, array_size, iblkx, iblkz;

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//


    // store an array of pef of the same size as nwinz
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    // loop over the blocks
    for (unsigned int iblk=0; iblk<blocks.size(); iblk++){


        xstart = blocks[iblk].xstart;
        zstart = blocks[iblk].zstart;
        xend = blocks[iblk].xend;
        zend = blocks[iblk].zend;
        iblkz = iblk - (iblk/nwinz) * nwinz;
        iblkx = iblk/nwinz;

        array_size = pef_array.size();

        // compute data norm2 in the current block
        norm2_dat = 0;
        for (int ix_d=xstart+nx_p-1; ix_d<xend; ix_d++){
                for (int iz_d=zstart+nz_p-1-_lead; iz_d<zend-_lead; iz_d++){
                    norm2_dat += (*_dat0->_mat)[ix_d][iz_d]*(*_dat0->_mat)[ix_d][iz_d];
                }
        }

        // pick the appropriate pef to start with
        if (iblk == 0)
            pef = pef_array[0];
        
        else if (iblkz == 0)
            pef = pef_array[1];

        else if (iblkx ==0)
            pef = pef_array[array_size-1];

        else{
            norm2_res = 0;
            norm2_res2 = 0;
            norm2_res3 = 0;
            norm2_res4 = 0;
            res_temp.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            res_temp2.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            res_temp3.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            res_temp4.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));

            for (int ix_r=nx_p-1; ix_r<xend-xstart; ix_r++){
                for (int iz_r=nz_p-1; iz_r<zend-zstart; iz_r++){
                    for (int ix_p=0; ix_p<nx_p; ix_p++){
                        for (int iz_p=0; iz_p<nz_p; iz_p++){
                            res_temp[ix_r][iz_r] += (*pef_array[array_size-1]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                            res_temp2[ix_r][iz_r] += (*pef_array[0]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                            res_temp3[ix_r][iz_r] += (*pef_array[1]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                            res_temp4[ix_r][iz_r] += (*pef_array[2]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                        }
                    }
                    norm2_res += res_temp[ix_r][iz_r]*res_temp[ix_r][iz_r];
                    norm2_res2 += res_temp2[ix_r][iz_r]*res_temp2[ix_r][iz_r];
                    norm2_res3 += res_temp3[ix_r][iz_r]*res_temp3[ix_r][iz_r];
                    norm2_res4 += res_temp4[ix_r][iz_r]*res_temp4[ix_r][iz_r];
                }
            }

            if (norm2_res >= norm2_dat) {w_up = 0;}
            else {w_up = 1 - norm2_res / norm2_dat;}

            if (norm2_res2 >= norm2_dat) {w_up_left = 0;}
            else {w_up_left = 1 - norm2_res2/ norm2_dat;}

            if (norm2_res3 >= norm2_dat) {w_left = 0;}
            else {w_left = 1 - norm2_res3/ norm2_dat;}

            if (norm2_res4 >= norm2_dat) {w_dn_left = 0;}
            else {w_dn_left = 1 - norm2_res4/ norm2_dat;}

            w_sum = w_up + w_up_left + w_left + w_dn_left;

            if (w_sum != 0){
                w_up /= w_sum;
                w_up_left /= w_sum;
                w_left /= w_sum;
                w_dn_left /= w_sum;
            }

            pef = pef_array[array_size-1];
            pef->scaleAdd(pef_array[0], w_up, w_up_left);
            pef->scaleAdd(pef_array[1], 1.0, w_left);
            pef->scaleAdd(pef_array[2], 1.0, w_dn_left);
        }

        (*pef->_mat)[0][_lead] = 1;

        // inner loop to estimate the pef using steepest descent
        for (int iter = 0; iter < _niter; iter++){

            // compute residuals
            res_temp.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int ix_d=xstart+nx_p-1-ix_p; ix_d<xend-ix_p; ix_d++){
                    for (int iz_p=0; iz_p<nz_p; iz_p++){
                        for (int iz_d=zstart+nz_p-1-iz_p; iz_d<zend-iz_p; iz_d++){
                            res_temp[ix_d-xstart+ix_p][iz_d-zstart+iz_p] += (*pef->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_d][iz_d];
                        }
                    }
                }
            }

            // compute filter perturbation
            dpef->zero();
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int ix_d=xstart+nx_p-1-ix_p; ix_d<xend-ix_p; ix_d++){
                    for (int iz_p=0; iz_p<nz_p; iz_p++){
                        for (int iz_d=zstart+nz_p-1-iz_p; iz_d<zend-iz_p; iz_d++){
                            (*dpef->_mat)[ix_p][iz_p] += res_temp[ix_d-xstart+ix_p][iz_d-zstart+iz_p]*(*_dat0->_mat)[ix_d][iz_d];
                        }
                    }
                }
            }

            // compute steepest descent step length
            norm2_dg = 0;
            res_temp.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            for (int ix_r=nx_p-1; ix_r<xend-xstart; ix_r++){
                for (int iz_r=nz_p-1; iz_r<zend-zstart; iz_r++){
                    for (int ix_p=0; ix_p<nx_p; ix_p++){
                        for (int iz_p=0; iz_p<nz_p; iz_p++){
                            res_temp[ix_r][iz_r] += (*dpef->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                        }
                    }
                    norm2_dg += res_temp[ix_r][iz_r]*res_temp[ix_r][iz_r];
                }
            }
            
            norm2_g = VECEXT::norm2(dpef);
            if ((norm2_g == 0) || (norm2_dg<1e-07) ){
                break;
            }
            else{
                step_length = norm2_g / norm2_dg;
            }

            // update the pef
            for (int iz_p=0; iz_p<=_lead; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
            }
            pef->scaleAdd(dpef, 1.0, - step_length);
        }

        // store the new pef
        newpef = pef->clone();
        pef_array.push_back(newpef);
        array_size = pef_array.size();

        // apply the new PEF to the output using 4 weighted PEFs, the current, the upper, the left and the upper left one
        norm2_res = 0;
        for (int ix_r=xstart+nx_p-1; ix_r<xend; ix_r++){
            for (int iz_r = zstart+nz_p-1-_lead; iz_r<zend-_lead; iz_r++){

                dx = (ix_r -xstart + 1)/_winx;
                dx = std::min(1, iblkx)*dx - std::min(0, iblkx - 1);
                dz = (iz_r - zstart - nz_p + _lead + 2)/_winz;
                dz = std::min(1, iblkz)*dz - std::min(0, iblkz - 1);

                assert (abs(dx)<=1);
                assert(abs(dz)<=1);
                assert(dx*dz+(dx-dx*dz)+(dz-dx*dz)+(1-dx)*(1-dz) == 1);

                for (int ix_p=0; ix_p<nx_p; ix_p++){
                    for (int iz_p=0; iz_p<nz_p; iz_p++){
                        (*dat->_mat)[ix_r][iz_r] += (*newpef->_mat)[ix_p][iz_p]* //current pef without interpolation
                                                    //(dx*dz*(*newpef->_mat)[ix_p][iz_p]+ // current pef
                                                    //(dx-dx*dz)*(*pef_array[array_size-2]->_mat)[ix_p][iz_p]+ // upper pef
                                                    //(dz-dx*dz)*(*pef_array[1]->_mat)[ix_p][iz_p]+ // left pef
                                                    //(1-dx)*(1-dz)*(*pef_array[0]->_mat)[ix_p][iz_p]) * //upper left pef
                                                    (*mod->_mat)[ix_r-ix_p][iz_r-iz_p+_lead];
                    }
                }
                norm2_res += (*dat->_mat)[ix_r][iz_r] * (*dat->_mat)[ix_r][iz_r];
            }
        }

        // remove the unwanted pef        
        if (iblk >= nwinz)
           pef_array.erase(pef_array.begin()); 
    }

    _peff = pef_array[pef_array.size()-1]->clone();
}



void nsPef2D_v9::adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const {
    this->checkDomainRange(mod,dat);

    if (add == false) mod->scale(0.0);

    // define arrays size
    int nx_p = _pef->getHyper()->getAxis(2).n;
    int nz_p = _pef->getHyper()->getAxis(1).n;
    int nx_d = mod->getHyper()->getAxis(2).n;
    int nz_d = mod->getHyper()->getAxis(1).n;

    // define a new structure for the local blocks
    struct block {
        int xstart, xend, zstart, zend, xsize, zsize;
    };

    // define the number of local blocks
    int nwinx = (nx_d - _winx) / (_winx - nx_p + 1);
    nwinx++;
    int nwinz = (nz_d - _winz) / (_winz - nz_p + 1);
    nwinz++;

//std::clog << "Number of blocks in x = " << nwinx << "\n";
//std::clog << "Number of blocks in z = " << nwinz << "\n";

    // build the list of local blocks
    std::vector<block> blocks (nwinx*nwinz);
    for (int ix=0; ix<nwinx; ix++){
        for (int iz=0; iz<nwinz; iz++){
            blocks[ix*nwinz+iz].xstart = ix*(_winx - nx_p + 1);
            blocks[ix*nwinz+iz].zstart = iz*(_winz - nz_p + 1);
            blocks[ix*nwinz+iz].xend = blocks[ix*nwinz+iz].xstart + _winx; 
            blocks[ix*nwinz+iz].zend = blocks[ix*nwinz+iz].zstart + _winz;

            if (ix == nwinx - 1){
                blocks[ix*nwinz+iz].xend = nx_d;
                }
            if (iz == nwinz - 1){
                blocks[ix*nwinz+iz].zend = nz_d;
            }

            blocks[ix*nwinz+iz].xsize = blocks[ix*nwinz+iz].xend - blocks[ix*nwinz+iz].xstart;
            blocks[ix*nwinz+iz].zsize = blocks[ix*nwinz+iz].zend - blocks[ix*nwinz+iz].zstart;
        }
    }

//std::clog << "Maximum block size in x = " << blocks[nwinx*nwinz-1].xsize << "\n";
//std::clog << "Maximum block size in z = " << blocks[nwinx*nwinz-1].zsize << "\n";

    // define residual vectors for the largest block (last one)
    std::vector<std::vector <float> > res_temp (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));
    std::vector<std::vector <float> > res_temp2 (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));
    std::vector<std::vector <float> > res_temp3 (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));
    std::vector<std::vector <float> > res_temp4 (blocks[nwinx*nwinz-1].xsize, std::vector<float>(blocks[nwinx*nwinz-1].zsize, 0));

    // define intermediate vectors and variables
    std::shared_ptr<float2DReg> pef;
    std::shared_ptr<float2DReg> newpef;
    std::shared_ptr<float2DReg> dpef = _pef->clone();
    dpef->zero();

    double step_length, norm2_dg, norm2_g, norm2_res, norm2_res2, norm2_res3, norm2_res4, norm2_dat, w_up, w_left, w_up_left, w_dn_left, w_sum, dx, dz;
    int xstart, xend, zstart, zend, array_size, iblkx, iblkz;

    //####################################################//
    // Estimation and application excluding the boundaries//
    //####################################################//


    // store an array of pef of the same size as nwinz
    std::vector<std::shared_ptr<float2DReg> > pef_array;
    pef_array.push_back(_pef->clone());

    // loop over the blocks
    for (unsigned int iblk=0; iblk<blocks.size(); iblk++){


        xstart = blocks[iblk].xstart;
        zstart = blocks[iblk].zstart;
        xend = blocks[iblk].xend;
        zend = blocks[iblk].zend;
        iblkz = iblk - (iblk/nwinz) * nwinz;
        iblkx = iblk/nwinz;

        array_size = pef_array.size();

        // compute data norm2 in the current block
        norm2_dat = 0;
        for (int ix_d=xstart+nx_p-1; ix_d<xend; ix_d++){
                for (int iz_d=zstart+nz_p-1-_lead; iz_d<zend-_lead; iz_d++){
                    norm2_dat += (*_dat0->_mat)[ix_d][iz_d]*(*_dat0->_mat)[ix_d][iz_d];
                }
        }


        // pick the appropriate pef to start with
        if (iblk == 0)
            pef = pef_array[0];
        
        else if (iblkz == 0)
            pef = pef_array[1];

        else if (iblkx ==0)
            pef = pef_array[array_size-1];

        else{
            norm2_res = 0;
            norm2_res2 = 0;
            norm2_res3 = 0;
            norm2_res4 = 0;
            res_temp.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            res_temp2.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            res_temp3.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            res_temp4.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));

            for (int ix_r=nx_p-1; ix_r<xend-xstart; ix_r++){
                for (int iz_r=nz_p-1; iz_r<zend-zstart; iz_r++){
                    for (int ix_p=0; ix_p<nx_p; ix_p++){
                        for (int iz_p=0; iz_p<nz_p; iz_p++){
                            res_temp[ix_r][iz_r] += (*pef_array[array_size-1]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                            res_temp2[ix_r][iz_r] += (*pef_array[0]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                            res_temp3[ix_r][iz_r] += (*pef_array[1]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                            res_temp4[ix_r][iz_r] += (*pef_array[2]->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                        }
                    }
                    norm2_res += res_temp[ix_r][iz_r]*res_temp[ix_r][iz_r];
                    norm2_res2 += res_temp2[ix_r][iz_r]*res_temp2[ix_r][iz_r];
                    norm2_res3 += res_temp3[ix_r][iz_r]*res_temp3[ix_r][iz_r];
                    norm2_res4 += res_temp4[ix_r][iz_r]*res_temp4[ix_r][iz_r];
                }
            }

            if (norm2_res >= norm2_dat) {w_up = 0;}
            else {w_up = 1 - norm2_res / norm2_dat;}

            if (norm2_res2 >= norm2_dat) {w_up_left = 0;}
            else {w_up_left = 1 - norm2_res2/ norm2_dat;}

            if (norm2_res3 >= norm2_dat) {w_left = 0;}
            else {w_left = 1 - norm2_res3/ norm2_dat;}

            if (norm2_res4 >= norm2_dat) {w_dn_left = 0;}
            else {w_dn_left = 1 - norm2_res4/ norm2_dat;}

            w_sum = w_up + w_up_left + w_left + w_dn_left;

            if (w_sum != 0){
                w_up /= w_sum;
                w_up_left /= w_sum;
                w_left /= w_sum;
                w_dn_left /= w_sum;
            }

            pef = pef_array[array_size-1];
            pef->scaleAdd(pef_array[0], w_up, w_up_left);
            pef->scaleAdd(pef_array[1], 1.0, w_left);
            pef->scaleAdd(pef_array[2], 1.0, w_dn_left);
        }

        (*pef->_mat)[0][_lead] = 1;

        // inner loop to estimate the pef using steepest descent
        for (int iter = 0; iter < _niter; iter++){

            // compute residuals
            res_temp.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int ix_d=xstart+nx_p-1-ix_p; ix_d<xend-ix_p; ix_d++){
                    for (int iz_p=0; iz_p<nz_p; iz_p++){
                        for (int iz_d=zstart+nz_p-1-iz_p; iz_d<zend-iz_p; iz_d++){
                            res_temp[ix_d-xstart+ix_p][iz_d-zstart+iz_p] += (*pef->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_d][iz_d];
                        }
                    }
                }
            }

            // compute filter perturbation
            dpef->zero();
            for (int ix_p=0; ix_p<nx_p; ix_p++){
                for (int ix_d=xstart+nx_p-1-ix_p; ix_d<xend-ix_p; ix_d++){
                    for (int iz_p=0; iz_p<nz_p; iz_p++){
                        for (int iz_d=zstart+nz_p-1-iz_p; iz_d<zend-iz_p; iz_d++){
                            (*dpef->_mat)[ix_p][iz_p] += res_temp[ix_d-xstart+ix_p][iz_d-zstart+iz_p]*(*_dat0->_mat)[ix_d][iz_d];
                        }
                    }
                }
            }

            // compute steepest descent step length
            norm2_dg = 0;
            res_temp.assign(res_temp.size(), std::vector<float>(res_temp[0].size(),0));
            for (int ix_r=nx_p-1; ix_r<xend-xstart; ix_r++){
                for (int iz_r=nz_p-1; iz_r<zend-zstart; iz_r++){
                    for (int ix_p=0; ix_p<nx_p; ix_p++){
                        for (int iz_p=0; iz_p<nz_p; iz_p++){
                            res_temp[ix_r][iz_r] += (*dpef->_mat)[ix_p][iz_p]*(*_dat0->_mat)[ix_r+xstart-ix_p][iz_r+zstart-iz_p];
                        }
                    }
                    norm2_dg += res_temp[ix_r][iz_r]*res_temp[ix_r][iz_r];
                }
            }
            
            norm2_g = VECEXT::norm2(dpef);
            if ((norm2_g == 0) || (norm2_dg<1e-07) ){
                break;
            }
            else{
                step_length = norm2_g / norm2_dg;
            }

            // update the pef
            for (int iz_p=0; iz_p<=_lead; iz_p++){
                (*dpef->_mat)[0][iz_p] = 0;
            }
            pef->scaleAdd(dpef, 1.0, - step_length);
        }

        // store the new pef
        newpef = pef->clone();
        pef_array.push_back(newpef);
        array_size = pef_array.size();

        // apply the new PEF to the output using 4 weighted PEFs, the current, the upper, the left and the upper left one
        norm2_res = 0;
        for (int ix_r=xstart+nx_p-1; ix_r<xend; ix_r++){
            for (int iz_r = zstart+nz_p-1-_lead; iz_r<zend-_lead; iz_r++){

                dx = (ix_r -xstart + 1)/_winx;
                dx = std::min(1, iblkx)*dx - std::min(0, iblkx - 1);
                dz = (iz_r - zstart - nz_p + _lead + 2)/_winz;
                dz = std::min(1, iblkz)*dz - std::min(0, iblkz - 1);

                assert (abs(dx)<=1);
                assert(abs(dz)<=1);
                assert(dx*dz+(dx-dx*dz)+(dz-dx*dz)+(1-dx)*(1-dz) == 1);

                for (int ix_p=0; ix_p<nx_p; ix_p++){
                    for (int iz_p=0; iz_p<nz_p; iz_p++){
                        (*mod->_mat)[ix_r-ix_p][iz_r-iz_p+_lead] += (*newpef->_mat)[ix_p][iz_p]* //current pef without interpolation
                                                    //(dx*dz*(*newpef->_mat)[ix_p][iz_p]+ // current pef
                                                    //(dx-dx*dz)*(*pef_array[array_size-2]->_mat)[ix_p][iz_p]+ // upper pef
                                                    //(dz-dx*dz)*(*pef_array[1]->_mat)[ix_p][iz_p]+ // left pef
                                                    //(1-dx)*(1-dz)*(*pef_array[0]->_mat)[ix_p][iz_p]) * //upper left pef
                                                    (*dat->_mat)[ix_r][iz_r];
                    }
                }
                norm2_res += (*mod->_mat)[ix_r][iz_r] * (*mod->_mat)[ix_r][iz_r];
            }
        }

        // remove the unwanted pef        
        if (iblk >= nwinz)
           pef_array.erase(pef_array.begin()); 
    }

    _peff = pef_array[pef_array.size()-1]->clone();
}
