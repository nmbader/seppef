#ifndef doubleHyperExt_H
#define doubleHyperExt_H

#include <iostream>
#include <string>
#include "double1DReg.h"
#include "double2DReg.h"
#include "doubleHyper.h"

using namespace SEP;
namespace VECDEXT {

// Non-member functions to extend functionalities of doubleHyper vectors 

// write to SEP format (.H file)
void sepWrite(const std::shared_ptr<doubleHyper> vec, std::string output="out");

// read from SEP format (.H file) given in command line by '<'
std::shared_ptr<double1DReg> sepRead1D();

// read from SEP format (.H file) given in command line by 'input=file.H'
std::shared_ptr<double1DReg> sepRead1D(std::string input);

// read from SEP format (.H file) given in command line by '<'
std::shared_ptr<double2DReg> sepRead();

// read from SEP format (.H file) given in command line by 'input=file.H'
std::shared_ptr<double2DReg> sepRead(std::string input);

// read from SEP format (.H file) given in command line by '<'
std::shared_ptr<double3DReg> sepRead3D();

// read from SEP format (.H file) given in command line by 'input=file.H'
std::shared_ptr<double3DReg> sepRead3D(std::string input);

// fill the vector with random numbers
void random(std::shared_ptr<doubleHyper> vec, float min=-1, float max=1);

// compute norm2 squared of the vector
double norm2(const std::shared_ptr<doubleHyper> vec);

// compute the sum of all samples
double sum(const std::shared_ptr<doubleHyper> vec);

// compute the RMS value of the vector
double rms(const std::shared_ptr<doubleHyper> vec);

} // VECDEXT

#endif