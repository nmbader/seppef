#include <time.h>
#include <chrono>

#include <sys/stat.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "axis.h"
#include "floatHyperExt.h"
#include "float2DReg.h"
#include "functions.h"
#include "convolution2D.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program convolve2D \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");
    

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in the 2D filter
    std::string fil_file;
    MB::readParam(argc, argv, "filter", fil_file,"none");
    if (fil_file == "none") throw std::logic_error("The required filter is not provided.\n");
    std::shared_ptr<float2DReg> fil = VECEXT::sepRead(fil_file);
    unsigned int fil_nx = fil->getHyper()->getAxis(2).n;
    unsigned int fil_nz = fil->getHyper()->getAxis(1).n;

// read in other parameters
    bool crossBounds, reshape, addBounds;
    MB::readParam(argc, argv, "cross_bounds", crossBounds, false);
    MB::readParam(argc, argv, "add_bounds", addBounds, false);
    MB::readParam(argc, argv, "reshape", reshape, true);

// convolve the filter with the data
    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);
    X.n = VECEXT::getnx(dat) + VECEXT::getnx(fil) - 1;
	Z.n = VECEXT::getnz(dat) + VECEXT::getnz(fil) - 1;
	X.o = VECEXT::getox(dat) + VECEXT::getox(fil);
	Z.o = VECEXT::getoz(dat) + VECEXT::getoz(fil);
    std::shared_ptr<float2DReg> dat0 (new float2DReg(Z,X));

    if ((crossBounds == true) || ((crossBounds == false) && (addBounds==false))){
        convolution2D conv(fil,dat->getHyper(),dat0->getHyper());
        conv.crossBounds(crossBounds);
        conv.forward(false, dat, dat0);
    }

    else{
        (*fil->_mat)[0][0] -= 1;
        convolution2D conv(fil,dat->getHyper(),dat0->getHyper());
        conv.crossBounds(crossBounds);
        conv.forward(false, dat, dat0);
        fil = VECEXT::resize(fil,1,1);
	    (*fil->_mat)[0][0] = 1;
        convolution2D conv0(fil,dat->getHyper(),dat0->getHyper());
        conv0.crossBounds(true);
        conv0.forward(true, dat, dat0); 
    }

//  reshape to the size of the input in case needed
    if (reshape == true){
        dat0 = VECEXT::reshape(dat0, dat);
    }

// write the result to disk
    VECEXT::sepWrite(dat0, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program convolve2D; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}