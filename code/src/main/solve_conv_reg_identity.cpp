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
#include "cglsReg.h"
#include "convolution2D.h"
#include "identity2D.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================

// solving min ||L.m - d|| + ||eps.(m - m_prior)||
// L is a filter convolution matrix
// m0 = input data


int main(int argc, char **argv) {
    
    std::clog << "Starting program: solving min|| L.m - d ||^2 + eps^2|| m - m_prior ||^2 ;\n";
    std::clog << "L = convolution with a filter \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");
    

// read in the input model
    std::shared_ptr<float2DReg> mod = VECEXT::sepRead();

// read in the filter
    std::string fil_file;
    MB::readParam(argc, argv, "fil", fil_file, "none");
    if (fil_file == "none") throw std::logic_error("The required filter is not provided.\n");
    std::shared_ptr<float2DReg> fil = VECEXT::sepRead(fil_file);
    unsigned int fil_nx = fil->getHyper()->getAxis(2).n;
    unsigned int fil_nz = fil->getHyper()->getAxis(1).n;

// read in other parameters
    int niter;
    float threshold, eps;
    bool cross_bounds;
    std::string objFunc, grad, mod_grad, dat_grad, mod_res, dat_res, solver;
    MB::readParam(argc, argv, "niter", niter, 100);
    MB::readParam(argc, argv, "threshold", threshold, 0.01);
    MB::readParam(argc, argv, "obj_func", objFunc,"none");
    MB::readParam(argc, argv, "grad", grad,"none");
    MB::readParam(argc, argv, "mod_grad", mod_grad,"none");
    MB::readParam(argc, argv, "dat_grad", dat_grad,"none");
    MB::readParam(argc, argv, "mod_res", mod_res,"none");
    MB::readParam(argc, argv, "dat_res", dat_res,"none");
    MB::readParam(argc, argv, "solver", solver,"cgls");
    MB::readParam(argc, argv, "eps", eps, 0.0);
    MB::readParam(argc, argv, "cross_bounds", cross_bounds,false);

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

    std::shared_ptr<float2DReg> prior;
    std::string prior_file;
    MB::readParam(argc, argv, "prior", prior_file,"none");
 
    if (prior_file != "none"){
        prior = VECEXT::sepRead(prior_file);
    }
    else{
        prior = nullptr;
    }
    
    convolution2D L(fil,mod->getHyper(),dat->getHyper()); // define the filter convolution matrix L
    L.crossBounds(cross_bounds);

    identity2D I(mod->getHyper()); // define the filter convolution operator D

// set the solver
    std::vector<std::string> list = {"cgls"}; // list of available solvers
    std::shared_ptr<lsolverReg> sol;
    if (solver == "cgls"){
        std::shared_ptr<cglsReg> cg (new cglsReg(niter,threshold));
        sol = cg;
    }

    else {
        throw std::logic_error("The requested solver is not implemented.\n");
    }

// solve the problem
    sol->run(&L, &I, eps, mod, dat, prior);

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

    if (mod_grad != "none"){
        VECEXT::sepWrite(sol->_m_grad, mod_grad);
    }

    if (dat_grad != "none"){
        std::shared_ptr<floatHyper> d_grad (new float2DReg(sol->_grad->getHyper()));
        d_grad->scaleAdd(sol->_grad,0,1);
        d_grad->scaleAdd(sol->_m_grad,1,-1);
        VECEXT::sepWrite(d_grad, dat_grad);
    }

    if (mod_res != "none"){
        VECEXT::sepWrite(sol->_m_res, mod_res);
    }

    if (dat_res != "none"){
        VECEXT::sepWrite(sol->_d_res, dat_res);
    }

// write the model to disk
    VECEXT::sepWrite(mod, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}