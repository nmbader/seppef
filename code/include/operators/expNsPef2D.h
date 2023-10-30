#ifndef expNsPef2D_H
#define expNsPef2D_H

#include "nsPef2D.h"

using namespace SEP;

// Non-stationary 2D PEF operators, expanded when applied


class expNsPef2D_v0 : public nsPef2D_v0{
public:
    // inherit all constructors from the base class
    using nsPef2D_v0::nsPef2D_v0;

    // destuctor
    ~expNsPef2D_v0(){}

    // forward and adjoint operators are the same as for the base class, except that the PEF
    // is expanded before application

    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
};

class expNsPef2D_v1 : public nsPef2D_v1{
public:
    // inherit all constructors from the base class
    using nsPef2D_v1::nsPef2D_v1;

    // destuctor
    ~expNsPef2D_v1(){}

    // forward and adjoint operators are the same as for the base class, except that the PEF
    // is expanded before application

    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
};

class expNsPef2D_v2 : public nsPef2D_v2{
public:
    // inherit all constructors from the base class
    using nsPef2D_v2::nsPef2D_v2;

    // destuctor
    ~expNsPef2D_v2(){}

    // forward and adjoint operators are the same as for the base class, except that the PEF
    // is expanded before application

    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
};

class expNsPef2D_v3 : public nsPef2D_v3{
public:
    // inherit all constructors from the base class
    using nsPef2D_v3::nsPef2D_v3;

    // destuctor
    ~expNsPef2D_v3(){}

    // forward and adjoint operators are the same as for the base class, except that the PEF
    // is expanded before application

    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
};

class expNsPef2D_v4 : public nsPef2D_v4{
public:
    // inherit all constructors from the base class
    using nsPef2D_v4::nsPef2D_v4;

    // destuctor
    ~expNsPef2D_v4(){}

    // forward and adjoint operators are the same as for the base class, except that the PEF
    // is expanded before application

    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
};

class expNsPef2D_v5 : public nsPef2D_v5{
public:
    // inherit all constructors from the base class
    using nsPef2D_v5::nsPef2D_v5;

    // destuctor
    ~expNsPef2D_v5(){}

    // forward and adjoint operators are the same as for the base class, except that the PEF
    // is expanded before application

    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
};

class expNsPef2D_v9 : public nsPef2D_v9{
public:
    // inherit all constructors from the base class
    using nsPef2D_v9::nsPef2D_v9;

    // destuctor
    ~expNsPef2D_v9(){}

    // forward and adjoint operators are the same as for the base class, except that the PEF
    // is expanded before application

    void forward(bool add, const std::shared_ptr<float2DReg> mod, std::shared_ptr<float2DReg> dat) const;
    void adjoint(bool add, std::shared_ptr<float2DReg> mod, const std::shared_ptr<float2DReg> dat) const;
};


#endif