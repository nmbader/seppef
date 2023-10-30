#ifndef chainOper2D_H
#define chainOper2D_H

#include "oper2D.h"
#include "float2DReg.h"

using namespace SEP;

// concatenation of two operators
class chainOper2D : public oper2D{
private:
    oper2D * _operLeft;
    oper2D * _operRight;

public:
    // default constructor
    chainOper2D(){_operLeft = nullptr; _operRight = nullptr;}

    // destructor
    ~chainOper2D(){delete _operLeft; delete _operRight;}

    // constructor from two operators
    chainOper2D(const oper2D* opLeft, const oper2D * opRight){
        opLeft->getDomain()->checkSame(opRight->getRange());
        _operLeft = opLeft->clone();
        _operRight = opRight->clone();
        _domain = opRight->getDomain()->clone();
        _range = opLeft->getRange()->clone();
    }

    // cloning the object
    chainOper2D * clone() const {
        chainOper2D * dOp = new chainOper2D(_operLeft, _operRight);
        return dOp;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const{
        std::shared_ptr<float2DReg> datRight (new float2DReg(_operRight->getRange()));
        _operRight->forward(add, mod, datRight);
        _operLeft->forward(add, datRight, dat);
    }

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const{
        std::shared_ptr<float2DReg> modLeft (new float2DReg(_operLeft->getDomain()));
        _operLeft->adjoint(add, modLeft, dat);
        _operRight->adjoint(add, mod, modLeft);
    }

};

#endif
