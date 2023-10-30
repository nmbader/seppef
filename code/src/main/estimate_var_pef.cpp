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
#include "varPef2D.h"
#include "cgls.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program estimate_var_pef \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in the 2D PEF parameters
    int pef_nx, pef_nz, pef_lead, xinc, zinc;
    MB::readParam(argc, argv, "n1", pef_nz, 6);
    MB::readParam(argc, argv, "n2", pef_nx, 3);
    MB::readParam(argc, argv, "lead", pef_lead,0);
    MB::readParam(argc, argv, "xinc", xinc,50);
    MB::readParam(argc, argv, "zinc", zinc,50);

// read in the CGLS parameters
int niter;
    float threshold;
    MB::readParam(argc, argv, "niter", niter, 9999);
    MB::readParam(argc, argv, "threshold", threshold, 0.0001);
    cgls cg(niter, threshold);

// read in the initial PEFs if provided
    std::string allpef_file;
	MB::readParam(argc, argv, "allpef", allpef_file, "none");
    std::shared_ptr<float2DReg> allpef0;
	if (allpef_file != "none"){
    	allpef0 = VECEXT::sepRead(allpef_file);
	}
    else {
        allpef0 = nullptr;
    }

// initialize the variable PEF
    varPef2D p(dat->getHyper(), pef_nx, pef_nz, xinc, zinc, pef_lead, allpef0);

// read in the initial PEF if provided (this will override allpef0)
    std::string pef_file;
	MB::readParam(argc, argv, "pef", pef_file, "none");
	if (pef_file != "none"){
    	std::shared_ptr<float2DReg> pef0 = VECEXT::sepRead(pef_file);
        p.populate(pef0);
	}

// estimate the pef
    p.estimate(dat, cg);
    
// write the normalized objective function, gradient and residuals to disk if requested
    std::string objFunc, grad, res;
    MB::readParam(argc, argv, "obj_func", objFunc,"none");
    MB::readParam(argc, argv, "grad", grad,"none");
    MB::readParam(argc, argv, "res", res,"none");

    if (objFunc != "none"){
        std::shared_ptr<float1DReg> func (new float1DReg(cg._func.size()));
        for (int i=0; i<cg._func.size(); i++){
            func->getVals()[i] = cg._func[i];
        }

        VECEXT::sepWrite(func, objFunc);
    }

    if (grad != "none"){
        VECEXT::sepWrite(cg._grad, grad);
    }

    if (res != "none"){
        VECEXT::sepWrite(cg._res, res);
    }

// write the pef to disk
    VECEXT::sepWrite(p._allpef, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program estimate_var_pef; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}