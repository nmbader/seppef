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
#include "expVarConv2D.h"
#include "seplib.h"
#include "sdls.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================

// solving min ||F.m - d||
// F is an expanded filter convolution matrix containing bank of filters
// m0 = input data


int main(int argc, char **argv) {
    
    std::clog << "Starting program: solving min|| F.m - d ||_2 ; F = expanded convolution operator with bank of filters \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");
    

// read in the input model
    std::shared_ptr<float2DReg> mod = VECEXT::sepRead();

// read in the filter parameters
    int nx, nz, lead, xinc, zinc;
    std::string fil_file;
    MB::readParam(argc, argv, "n2", nx, 1);
    MB::readParam(argc, argv, "n1", nz, 1);
    MB::readParam(argc, argv, "lead", lead,0);
    MB::readParam(argc, argv, "xinc", xinc,1);
    MB::readParam(argc, argv, "zinc", zinc,1);
    MB::readParam(argc, argv, "fil", fil_file,"none");
    if (fil_file == "none") throw std::logic_error("The required filter is not provided.\n");
    std::shared_ptr<float2DReg> fil = VECEXT::sepRead(fil_file);

// read in other parameters
    int niter;
    float threshold;
    std::string objFunc, grad, res, solver;
    MB::readParam(argc, argv, "niter", niter, 100);
    MB::readParam(argc, argv, "threshold", threshold, 0.01);
    MB::readParam(argc, argv, "obj_func", objFunc,"none");
    MB::readParam(argc, argv, "grad", grad,"none");
    MB::readParam(argc, argv, "res", res,"none");
    MB::readParam(argc, argv, "solver", solver,"cgls");

// set the data vector
    std::shared_ptr<float2DReg> zero (new float2DReg(mod->getHyper())); // vector filled with zeros
    zero->zero();

    std::shared_ptr<float2DReg> dat;
    std::string dat_file;
    MB::readParam(argc, argv, "dat", dat_file,"none");
    if (dat_file != "none"){
        dat = VECEXT::sepRead(dat_file);
        dat = VECEXT::reshape(dat, zero);
    }
    else{
        dat = zero;
    }
    
    expVarConv2D conv(fil, mod->getHyper(), dat->getHyper(),
    nx, nz, xinc, zinc, lead); // define the filter convolution operator

// set the solver
    std::vector<std::string> list = {"cgls", "sdls"}; // list of available solvers
    std::shared_ptr<lsolver> sol;
    if (solver == "cgls"){
        std::shared_ptr<cgls> cg (new cgls(niter,threshold));
        sol = cg;
    }
    else if (solver=="sdls"){
        std::shared_ptr<sdls> sd (new sdls(niter,threshold));
        sol = sd;
    }
    else {
        throw std::logic_error("The requested solver is not implemented.\n");
    }

// solve the problem
    sol->run(&conv, mod, dat);

// write the normalized objective function, gradient and residuals to disk if requested
    if (objFunc != "none"){
        std::shared_ptr<float1DReg> func (new float1DReg(sol->_func.size()));
        for (int i=0; i<sol->_func.size(); i++){
            func->getVals()[i] = sol->_func[i];
        }

        VECEXT::sepWrite(func, objFunc);
    }

    if (grad != "none"){
        VECEXT::sepWrite(sol->_grad, grad);
    }

    if (res != "none"){
        VECEXT::sepWrite(sol->_res, res);
    }

// write the model to disk
    VECEXT::sepWrite(mod, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}