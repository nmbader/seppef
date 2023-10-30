#ifndef floatHyperExt_H
#define floatHyperExt_H

#include <iostream>
#include <string>
#include "float1DReg.h"
#include "float2DReg.h"
#include "floatHyper.h"
#include "doubleHyper.h"
#include "double2DReg.h"

using namespace SEP;
namespace VECEXT {

// Non-member functions to extend functionalities of floatHyper vectors 

// get ox
float getox(const std::shared_ptr<float2DReg> vec);

// get oz
float getoz(const std::shared_ptr<float2DReg> vec);

// get dx
float getdx(const std::shared_ptr<float2DReg> vec);

// get dz
float getdz(const std::shared_ptr<float2DReg> vec);

// get nx
unsigned int getnx(const std::shared_ptr<floatHyper> vec);

// get nz
unsigned int getnz(const std::shared_ptr<floatHyper> vec);

// set ox
void setox(std::shared_ptr<floatHyper> vec, float ox);

// set oz
void setoz(std::shared_ptr<floatHyper> vec, float oz);

// write data into binary file, plus an ascii header file of the same name
void write(const std::shared_ptr<floatHyper> vec, std::string name);

// read hypercube from an ascii header file before reading the data from binary file
std::shared_ptr<hypercube> readHyper(std::string name);

// read data from binary file
void read(const std::shared_ptr<floatHyper> vec, std::string name);

// write to SEP format (.H file)
void sepWrite(const std::shared_ptr<floatHyper> vec, std::string output="out");
void sepWrite(const std::shared_ptr<doubleHyper> vec, std::string output="out");

// read from SEP format (.H file) given in command line by '<'
std::shared_ptr<float1DReg> sepRead1D();

// read from SEP format (.H file) given in command line by 'input=file.H'
std::shared_ptr<float1DReg> sepRead1D(std::string input);

// read from SEP format (.H file) given in command line by '<'
std::shared_ptr<float2DReg> sepRead();

// read from SEP format (.H file) given in command line by 'input=file.H'
std::shared_ptr<float2DReg> sepRead(std::string input);

// read from SEP format (.H file) given in command line by '<'
std::shared_ptr<float3DReg> sepRead3D();

// read from SEP format (.H file) given in command line by 'input=file.H'
std::shared_ptr<float3DReg> sepRead3D(std::string input);

// copy the hypercube from another hypercube without the size
void copyHyper(std::shared_ptr<float2DReg> vec, const std::shared_ptr<hypercube> hyper);

// add a constant to the vector
void add(std::shared_ptr<float2DReg> vec, const float alpha);

// set all entries to a given value
void set(std::shared_ptr<floatHyper> vec, const float val);
void set(std::shared_ptr<doubleHyper> vec, const double val);

// invert the vector, using a threshold eps to avoid dividing by zero
void inverse(std::shared_ptr<float2DReg> vec, const float eps = 1e-07);

// revert the order of the samples in one or two dimensions
void revert(std::shared_ptr<float2DReg> vec, const bool xRevert = true, const bool zRevert = true);
void revert(std::shared_ptr<double2DReg> vec, const bool xRevert = true, const bool zRevert = true);

// fill the vector with random numbers
void random(std::shared_ptr<floatHyper> vec, float min=-1, float max=1);
void random(std::shared_ptr<doubleHyper> vec, float min=-1, float max=1);

// compute norm2 squared of the vector
double norm2(const std::shared_ptr<floatHyper> vec);
double norm2(const std::shared_ptr<doubleHyper> vec);

// compute the sum of all samples
double sum(const std::shared_ptr<floatHyper> vec);
double sum(const std::shared_ptr<doubleHyper> vec);

// compute the RMS value of the vector
double rms(const std::shared_ptr<floatHyper> vec);
double rms(const std::shared_ptr<doubleHyper> vec);

// transpose the vector by exchanging x and z axes
std::shared_ptr<float2DReg> transpose(const std::shared_ptr<float2DReg> vec);

// expand the vector in 1 or 2 dimensions by adding zeros
std::shared_ptr<float2DReg> expand(const std::shared_ptr<float2DReg> vec, bool xExpand = true, bool zExpand = true);

// interpolate the vector in 1 or 2 dimensions (sinc interpolation in frequency domain)
std::shared_ptr<float2DReg> interpolate(const std::shared_ptr<float2DReg> vec, bool xInterp, bool zInterp);

// resamples the vector in 1 or 2 dimensions (by dropping every other sample)
std::shared_ptr<float2DReg> resample(const std::shared_ptr<float2DReg> vec, bool xResamp, bool zResamp);

// resize the vector to new x-z lengths
std::shared_ptr<float2DReg> resize(const std::shared_ptr<float2DReg> vec, unsigned int nx, unsigned int nz);

// pad with zeros at the start and at the end
std::shared_ptr<float1DReg> pad(const std::shared_ptr<float1DReg> vec, unsigned int nz_start, unsigned int nz_end);

// pad with zeros at the start and at the end
std::shared_ptr<float2DReg> pad(const std::shared_ptr<float2DReg> vec, unsigned int nx_start, unsigned int nx_end,
                                unsigned int nz_start, unsigned int nz_end);

// remove samples at the start and at the end
std::shared_ptr<float1DReg> cut(const std::shared_ptr<float1DReg> vec, unsigned int nz_start, unsigned int nz_end);

std::shared_ptr<float2DReg> cut(const std::shared_ptr<float2DReg> vec, unsigned int nx_start, unsigned int nx_end,
                                unsigned int nz_start, unsigned int nz_end);

// reshape the 2D vector to the size of another vector; where they overlap keep the first, where they don't, keep the second
std::shared_ptr<float2DReg> reshape(const std::shared_ptr<float2DReg> vec1, const std::shared_ptr<float2DReg> vec2);

// average the 2D vector along x direction
std::shared_ptr<float1DReg> stackx(const std::shared_ptr<float2DReg> vec);

// average the 2D vector along z direction
std::shared_ptr<float1DReg> stackz(const std::shared_ptr<float2DReg> vec); 

// Apply windowed sinc low pass filter
// wc is the cutoff angular frequency (percentage of PI)
// N is the half length of the filter, in number of samples
// alpha is the coefficient used for cosine windowing
void lpSinc(std::shared_ptr<float2DReg> vec, float wc, int halfLength=13, float alpha=0.5);

// Apply IIR Butterworth low or high pass filter
// wc is the cutoff angular frequency (percentage of PI)
// lowPass = false gives high pass filter
// halfOrder is the half order of the filter
// visit http://www.kwon3d.com/theory/filtering/butt.html for equations
// visit http://www.ee.ic.ac.uk/hp/staff/dmb/courses/DSPDF/00800_TransIIR.pdf
// for filters transformation (slide 66)
void iirButterworth(std::shared_ptr<float1DReg> vec, float wc, bool lowPass=true, unsigned int halfOrder=2);
void iirButterworth(std::shared_ptr<float2DReg> vec, float wc, bool lowPass=true, unsigned int halfOrder=2);

// Apply IIR Butterworth band-pass filter
// wcLow is the low cutoff angular frequency (percentage of Nyquist)
// wcHigh is the high cutoff angular frequency (percentage of Nyquist)
// halfOrder is the half order of the filter. The effective order is 4x that.
void iirBpButterworth(std::shared_ptr<float1DReg> vec, float wcLow, float wcHigh, unsigned int halfOrder=2);
void iirBpButterworth(std::shared_ptr<float2DReg> vec, float wcLow, float wcHigh, unsigned int halfOrder=2);

// Apply F-K filtering with cosine taper
// kLow is the low cutoff wavenumber (percentage of Nyquist)
// kHigh is the high cutoff wavenumber (percentage of Nyquist)
// taper is the taper size outside the cutoff (percentage of Nyquist)
void fkFilter(std::shared_ptr<float2DReg> vec, float kLow, float kHigh, float taper);

// Transform the 1D vector into a zero-phase equivalent, centered at time 0 with odd number of samples
std::shared_ptr<float1DReg> zero_phase(const std::shared_ptr<float1DReg> vec);

// Transform the 1D vector into a minimum-phase equivalent using Kolgomoroff factorization
std::shared_ptr<float1DReg> minimum_phase(const std::shared_ptr<float1DReg> vec, float eps=1e-07);

// Apply 1D alpha trim smoothing. 0 <= alpha <=1. alpha = 0 <=> mean smoothing. alpha = 1 <=> median smoothing.
// refer to http://www.librow.com/articles/article-7 for a detailed method with edges treatment
// the code below is far from optimal
std::shared_ptr<float1DReg> alphaTrim(const std::shared_ptr<float1DReg> vec, unsigned int halfLength, float alpha=0.0);

std::shared_ptr<float2DReg> alphaTrim(const std::shared_ptr<float2DReg> vec, unsigned int halfLengthx, unsigned int halfLengthz, float alpha=0.0);

} // VECEXT

#endif