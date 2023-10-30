#ifndef functions_H
#define functions_H

#include "float1DReg.h"
#include "float2DReg.h"
#include "complex2DReg.h"
#include <math.h>
#include <complex>
#include <string>

using namespace SEP;

// Functions declaration

// Overload the + operator between float and complex
inline std::complex<float> operator + (const float a, const std::complex<float> z){
    std::complex<float> c;
    c.real(a+z.real());
    c.imag(z.imag());
    return c;
}

// Overload the * operator between float and complex
inline std::complex<float> operator * (const float a, const std::complex<float> z){
    std::complex<float> c;
    c.real(a*z.real());
    c.imag(a*z.imag());
    return c;
}

// Overload the / operator between complex and float
inline std::complex<float> operator / (const std::complex<float> z, const float a){
    std::complex<float> c;
    c.real(z.real()/a);
    c.imag(z.imag()/a);
    return c;
}

// Define the 'sinc' function
inline double sinc(double x){
    if (x==0.0) return 1.0;
    else return sin(x)/x;
}


namespace MB {

// display the values of the 1D vector
void print(std::shared_ptr<float1DReg> vec, std::string message="", std::ostream &o=std::cerr);

// display the values of the 2D vector
void print(std::shared_ptr<float2DReg> vec, std::string message="", std::ostream &o=std::cerr);

// display the values of the 2D complex vector
void print(std::shared_ptr<complex2DReg> cvec, std::string message="", std::ostream &o=std::cerr);

// Build a cosine windowed sinc wavelet
// wc is the cutoff angular frequency (percentage of PI)
// N is the half length of the filter, in number of samples
// alpha is the coefficient used for cosine windowing
std::shared_ptr<float1DReg> sincWavelet(float wc, int halfLength=13, float alpha=0.5);

// Build a centered Gaussian wavelet
// sig is the std deviation in number of samples
// N is the half length of the wavelet, in number of samples
std::shared_ptr<float1DReg> gaussianWavelet(float sig=10, int N=100);

// Build a centered Ricker wavelet
// wc is the central angular frequency (percentage of PI)
// sig is the std deviation in number of samples
// N is the half length of the wavelet, in number of samples
// Central freq is wc = sqrt(2)/(sig*si) where si is the sampling interval
std::shared_ptr<float1DReg> rickerWavelet(float wc, int N=100);

// read parameter from a list of strings and store its value into a variable
// the parameter string must be of type "param=value"
// this function is used to read command line arguments
void readParam(int argc, char **argv, std::string param, int &value, int defaultVal);

void readParam(int argc, char **argv, std::string param, float &value, float defaultVal);

void readParam(int argc, char **argv, std::string param, double &value, double defaultVal);

void readParam(int argc, char **argv, std::string param, std::string &value, std::string defaultVal);

void readParam(int argc, char **argv, std::string param, bool &value, bool defaultVal);

void readParam(int argc, char **argv, std::string param, std::vector<float> &value, float defaultVal);

void readParam(int argc, char **argv, std::string param, std::vector<double> &value, double defaultVal);

} // MB

#endif