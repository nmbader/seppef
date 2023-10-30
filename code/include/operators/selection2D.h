#ifndef selection2D_H
#define selection2D_H

#include "oper2D.h"

using namespace SEP;

// selector operator (domain = range)
class selection2D : public oper2D{
private:
    std::shared_ptr<float2DReg> _selector;

public:
    // default constructor
    selection2D(){_selector = nullptr;}

    // destructor
    ~selection2D(){}

    // constructor
    selection2D(const std::shared_ptr<float2DReg> select){
        unsigned int nx = select->getHyper()->getAxis(2).n;
        unsigned int nz = select->getHyper()->getAxis(1).n;

        for (int ix=0; ix<nx; ix++){
            for (int iz=0; iz<nz; iz++){
                if (((*select->_mat)[ix][iz]!=0) && 
                    ((*select->_mat)[ix][iz]!=1)){
                        throw std::invalid_argument("The provided float2DReg is not a selector");
                    }
            }
        }
        
        _selector = select->clone();
        _domain = select->getHyper()->clone();
        _range = select->getHyper()->clone();
    }

    // cloning the object
    selection2D * clone() const{
        selection2D * sel = new selection2D(_selector);
        return sel;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;

    // reverse the selector
    void reverse();

};

#endif