#include <time.h>
#include <chrono>

#include <sys/stat.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "floatHyperExt.h"
#include "float1DReg.h"
#include "float2DReg.h"
#include "functions.h"
#include "pef2D.h"
#include "cgls.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program estimate_pef \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in the 2D PEF parameters
    int pef_nx, pef_nz, pef_lead, pef_gap, niter;
    std::string objFunc;
    MB::readParam(argc, argv, "n1", pef_nz, 6);
    MB::readParam(argc, argv, "n2", pef_nx, 3);
    MB::readParam(argc, argv, "lead", pef_lead,0);
    MB::readParam(argc, argv, "gap", pef_gap,0);
    MB::readParam(argc, argv, "obj_func", objFunc,"none");
    MB::readParam(argc, argv, "niter", niter, 9999);

// read in the CGLS threshold and condition of crossing bounds or not
    float threshold;
    bool crossBounds;
    MB::readParam(argc, argv, "threshold", threshold, 0.0001);
    MB::readParam(argc, argv, "cross_bounds", crossBounds, false);

// read in the data mask if provided
    std::string mask_file;
    std::shared_ptr<float2DReg> mask;
	MB::readParam(argc, argv, "mask", mask_file, "none");
	if (mask_file != "none"){
    	mask = VECEXT::sepRead(mask_file);
	}

// estimate the pef
	pef2D p(pef_lead, pef_gap);
	std::shared_ptr<float2DReg> pef (new float2DReg(pef_nz,pef_nx));

    cgls cg(niter, threshold);

    if (mask_file == "none"){
	    p.estimate(pef, dat, cg, crossBounds, niter);
    }
    else{
        p.estimate(pef, dat, mask, cg);
    }

// write the normalized objective function to disk if requested
    if (objFunc != "none"){
        std::shared_ptr<float1DReg> func (new float1DReg(cg._func.size()));
        for (int i=0; i<cg._func.size(); i++){
            func->getVals()[i] = cg._func[i];
        }

        VECEXT::sepWrite(func, objFunc);
    }

// write the pef to disk
    VECEXT::sepWrite(pef, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program estimate_pef; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}