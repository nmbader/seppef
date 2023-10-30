#ifndef complexHyperExt_H
#define complexHyperExt_H

#include <iostream>
#include <string>
#include <complex>
#include "float2DReg.h"
#include "complex2DReg.h"


using namespace SEP;
namespace CVECEXT {

// Non-member functions to extend functionalities of floatHyper vectors

// get nx
unsigned int getnx(const std::shared_ptr<complex2DReg> cvec);

// get nz
unsigned int getnz(const std::shared_ptr<complex2DReg> cvec);

// constructor of complex2DReg from a float2DReg (creates a pure real complex vector)
std::shared_ptr<complex2DReg> makeComplex(const std::shared_ptr<float2DReg> vec);

// add a real constant to the vector
void add(std::shared_ptr<complex2DReg> cvec, const float z);

// scale the vector by a real constant
void scale(std::shared_ptr<complex2DReg> cvec, const float z);

// add a complex constant to the vector
void add(std::shared_ptr<complex2DReg> cvec, const std::complex<float> z);

// scale the vector by a complex constant
void scale(std::shared_ptr<complex2DReg> cvec, const std::complex<float> z);

// add another complex vector to the vector
void add(std::shared_ptr<complex2DReg> cvec1, const std::shared_ptr<complex2DReg> cvec2);

// scale the vector by a complex constant then add another scaled vector to it
void scaleAdd(std::shared_ptr<complex2DReg> cvec1, const std::shared_ptr<complex2DReg> cvec2, const std::complex<float> alpha, const std::complex<float> beta);

// fill the vector with random numbers
void random(std::shared_ptr<complex2DReg> cvec, float min=-1, float max=1);

// perform dot product with another vector
std::complex<double> dot(const std::shared_ptr<complex2DReg> cvec1, const std::shared_ptr<complex2DReg> cvec2);

// compute norm2 squared of the vector
double norm2(const std::shared_ptr<complex2DReg> cvec);

// complex conjugate of the vector 
void conjugate(std::shared_ptr<complex2DReg> cvec);

// complex conjugate transpose (adjoint) of the vector by exchanging x and z axes
std::shared_ptr<complex2DReg> adjoint(const std::shared_ptr<complex2DReg> cvec);

// output real part into float2DReg
std::shared_ptr<float2DReg> real(const std::shared_ptr<complex2DReg> cvec);

// output imaginary part into float2DReg
std::shared_ptr<float2DReg> imag(const std::shared_ptr<complex2DReg> cvec);

// output complex module into float2DReg
std::shared_ptr<float2DReg> module(const std::shared_ptr<complex2DReg> cvec);

// output complex module squared into float2DReg
std::shared_ptr<float2DReg> module2(const std::shared_ptr<complex2DReg> cvec);

// output complex phase (argument) into float2DReg
std::shared_ptr<float2DReg> phase(const std::shared_ptr<complex2DReg> cvec);

// resize the vector to new x-z lengths
std::shared_ptr<complex2DReg> resize(const std::shared_ptr<complex2DReg> vec, int nx, int nz);

// pad with zeros at the start and at the end
std::shared_ptr<complex2DReg> pad(const std::shared_ptr<complex2DReg> vec, int nx_start, int nx_end,
                                unsigned int nz_start, unsigned int nz_end);

// copy the hypercube from another hypercube without the size
void copyHyper(std::shared_ptr<complex2DReg> vec, const std::shared_ptr<hypercube> hyper);

} // CVECEXT

#endif