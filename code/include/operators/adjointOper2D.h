#ifndef adjointOper2D_H
#define adjointOper2D_H

#include "oper2D.h"
#include "float2DReg.h"

using namespace SEP;

// the adjoint of an operator
class adjointOper2D : public oper2D{
private:
    oper2D * _oper;

public:
    // default constructor
    adjointOper2D(){_oper = nullptr;}

    // destructor
    ~adjointOper2D(){delete _oper;}

    // constructor from an operator
    adjointOper2D(const oper2D* op){
        _oper = op->clone();
        _domain = op->getRange()->clone();
        _range = op->getDomain()->clone();
    }

    // cloning the object
    adjointOper2D * clone() const {
        adjointOper2D * op = new adjointOper2D(_oper);
        return op;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const{
        _oper->adjoint(add, dat, mod);
    }

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const{
        _oper->forward(add, dat, mod);
    }

};

#endif
