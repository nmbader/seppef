#ifndef identity1D_H
#define identity1D_H

#include "oper1D.h"
#include "hypercube.h"

using namespace SEP;

// identity1D operator (domain = range)
class identity1D : public oper1D{
private:

public:
    // default constructor
    identity1D(){}

    // destructor
    ~identity1D(){}

    // constructor
    identity1D(const std::shared_ptr<hypercube> sh){
        _domain = sh->clone();
        _range = sh->clone();
    }

    // cloning the object
    identity1D * clone() const{
        identity1D * id = new identity1D(_domain);
        return id;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float1DReg> mod, std::shared_ptr<float1DReg> dat) const;

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float1DReg> mod, const std::shared_ptr<float1DReg> dat) const;
    
};

#endif