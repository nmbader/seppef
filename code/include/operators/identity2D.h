#ifndef identity2D_H
#define identity2D_H

#include "oper2D.h"
#include "hypercube.h"

using namespace SEP;

// identity2D operator (domain = range)
class identity2D : public oper2D{
private:

public:
    // default constructor
    identity2D(){}

    // destructor
    ~identity2D(){}

    // constructor
    identity2D(const std::shared_ptr<hypercube> sh){
        _domain = sh->clone();
        _range = sh->clone();
    }

    // cloning the object
    identity2D * clone() const{
        identity2D * id = new identity2D(_domain);
        return id;
    }

    // forward operator
    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;

    // adjoint operator
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
    
};

#endif