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
#include "pef2D.h"
#include "cgls.h"
#include "convolution2D.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================

// solving min ||Pm|| by CGLS
// P is a filter convolution matrix
// m0 = input data


int main(int argc, char **argv) {
    
    std::clog << "Starting program pef_lf_denoise_v1 \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");
    

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in the 2D PEF
    std::string pef_file;
    MB::readParam(argc, argv, "pef", pef_file);
    std::shared_ptr<float2DReg> pef = VECEXT::sepRead(pef_file);
    unsigned int pef_nx = pef->getHyper()->getAxis(2).n;
    unsigned int pef_nz = pef->getHyper()->getAxis(1).n;

// read in other parameters
    int niter;
    float threshold;
    bool crossBounds;
    std::string objFunc, grad, res;
    MB::readParam(argc, argv, "niter", niter, 100);
    MB::readParam(argc, argv, "threshold", threshold, 0.01);
    MB::readParam(argc, argv, "cross_bounds", crossBounds, false);
    MB::readParam(argc, argv, "obj_func", objFunc,"none");

// estimate the model (de-noised data)
    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);
    X.n = VECEXT::getnx(dat) + VECEXT::getnx(pef) - 1;
	Z.n = VECEXT::getnz(dat) + VECEXT::getnz(pef) - 1;
	X.o = VECEXT::getox(dat) + VECEXT::getox(pef);
	Z.o = VECEXT::getoz(dat) + VECEXT::getoz(pef);
    std::shared_ptr<float2DReg> zero (new float2DReg(Z,X)); // vector filled with zeros
    zero->zero();
    convolution2D conv(pef,dat->getHyper(),zero->getHyper()); // define the PEF convolution matrix
    conv.crossBounds(crossBounds);

    cgls cg(niter,threshold);
	cg.run(&conv, dat, zero);

// write the normalized objective function to disk if requested
    if (objFunc != "none"){
        std::shared_ptr<float1DReg> func (new float1DReg(cg._func.size()));
        for (int i=0; i<cg._func.size(); i++){
            func->getVals()[i] = cg._func[i];
        }

        VECEXT::sepWrite(func, objFunc);
    }

// write the model to disk
    VECEXT::sepWrite(dat, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program pef_lf_denoise_v1; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}