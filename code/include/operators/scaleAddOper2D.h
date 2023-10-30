#ifndef scaleAddOper2D_H
#define scaleAddOper2D_H

#include "oper2D.h"
#include "float2DReg.h"

using namespace SEP;

// scaled sum of two operators L + alpha*R
class scaleAddOper2D : public oper2D{
private:
    oper2D * _operLeft;
    oper2D * _operRight;
    float _scaleLeft;
    float _scaleRight;

public:
    // default constructor
    scaleAddOper2D(){_operLeft = nullptr; _operRight = nullptr;}

    // destructor
    ~scaleAddOper2D(){delete _operLeft; delete _operRight;}

    // constructor from two operators
    scaleAddOper2D(const oper2D* opLeft, const oper2D * opRight, const float scaleLeft = 1.0, const float scaleRight = 1.0){
        opLeft->checkDomainRange(opRight->getDomain(), opRight->getRange());
        assert(scaleLeft != 0);
        _operLeft = opLeft->clone();
        _operRight = opRight->clone();
        _domain = opRight->getDomain()->clone();
        _range = opLeft->getRange()->clone();
        _scaleLeft = scaleLeft;
        _scaleRight = scaleRight;
    }

    // cloning the object
    scaleAddOper2D * clone() const {
        scaleAddOper2D * dOp = new scaleAddOper2D(_operLeft, _operRight, _scaleLeft, _scaleRight);
        return dOp;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const{
        std::shared_ptr<float2DReg> temp (new float2DReg(dat->getHyper()));
        dat->scale(1.0/_scaleLeft);
        _operLeft->forward(add, mod, dat);
        _operRight->forward(false, mod, temp);
        dat->scaleAdd(temp, _scaleLeft, _scaleRight);
    }

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const{
        std::shared_ptr<float2DReg> temp (new float2DReg(mod->getHyper()));
        mod->scale(1.0/_scaleLeft);
        _operLeft->adjoint(add, mod, dat);
        _operRight->adjoint(false,temp,dat);
        mod->scaleAdd(temp, _scaleLeft, _scaleRight);
    }

};

#endif
