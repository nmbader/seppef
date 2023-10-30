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
#include "sdls.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================

// solving min ||F.m - d||
// F is a filter convolution matrix
// m0 = input data


int main(int argc, char **argv) {
    
    std::clog << "Starting program: solving min|| F.m - d ||_2 ; F = convolution operator \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");
    

// read in the input model
    std::shared_ptr<float2DReg> mod = VECEXT::sepRead();

// read in the filter
    std::string fil_file;
    MB::readParam(argc, argv, "fil", fil_file,"none");
    if (fil_file == "none") throw std::logic_error("The required filter is not provided.\n");
    std::shared_ptr<float2DReg> fil = VECEXT::sepRead(fil_file);
    unsigned int fil_nx = fil->getHyper()->getAxis(2).n;
    unsigned int fil_nz = fil->getHyper()->getAxis(1).n;

// read in other parameters
    int niter;
    float threshold;
    bool crossBounds;
    std::string objFunc, grad, res, solver;
    MB::readParam(argc, argv, "niter", niter, 100);
    MB::readParam(argc, argv, "threshold", threshold, 0.01);
    MB::readParam(argc, argv, "cross_bounds", crossBounds, false);
    MB::readParam(argc, argv, "obj_func", objFunc,"none");
    MB::readParam(argc, argv, "grad", grad,"none");
    MB::readParam(argc, argv, "res", res,"none");
    MB::readParam(argc, argv, "solver", solver,"cgls");

// set the data vector
    axis X = mod->getHyper()->getAxis(2);
    axis Z = mod->getHyper()->getAxis(1);
    X.n = VECEXT::getnx(mod) + VECEXT::getnx(fil) - 1;
	Z.n = VECEXT::getnz(mod) + VECEXT::getnz(fil) - 1;
	X.o = VECEXT::getox(mod) + VECEXT::getox(fil);
	Z.o = VECEXT::getoz(mod) + VECEXT::getoz(fil);
    std::shared_ptr<float2DReg> zero (new float2DReg(Z,X)); // vector filled with zeros
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
    
    convolution2D conv(fil,mod->getHyper(),dat->getHyper()); // define the filter convolution matrix
    conv.crossBounds(crossBounds);

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