#include <time.h>
#include <chrono>

#include <sys/stat.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "floatHyperExt.h"
#include "float2DReg.h"
#include "functions.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program fk_filter \n";
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();
    std::shared_ptr<float2DReg> temp (new float2DReg(dat->getHyper()));
    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);

// revert the X axis
    for (int ix=0; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*temp->_mat)[X.n-ix-1][iz] = (*dat->_mat)[ix][iz];
        }
    }

// read in the filter parameters
    float kLow, kHigh, taper;
    MB::readParam(argc, argv, "kmin", kLow, -0.5);
    MB::readParam(argc, argv, "kmax", kHigh, 0.5);
    MB::readParam(argc, argv, "taper", taper, 0.1);

// filter the data
    VECEXT::fkFilter(temp, kLow, kHigh, taper);

// revert back the X axis
    for (int ix=0; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*dat->_mat)[X.n-ix-1][iz] = (*temp->_mat)[ix][iz];
        }
    }   

// write the filtered data to disk
    VECEXT::sepWrite(dat, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program fk_filter; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}