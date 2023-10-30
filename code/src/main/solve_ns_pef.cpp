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
#include "nsPef2D.h"
#include "lsolver.h"
#include "cgls.h"
#include "sdls.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================

// solving min ||Ns.m - d||
// Ns is a non-stationary PEF
// m0 = input data


int main(int argc, char **argv) {
    
    std::clog << "Starting program: solving min|| Ns.m - d ||_2 ; Ns = non-stationary PEF \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data from which to estimate the PEF
    std::shared_ptr<float2DReg> dat0 = VECEXT::sepRead();

// read in the data to which the PEF is to be applied
    std::string mod_file;
    MB::readParam(argc, argv, "input", mod_file,"none");
    if (mod_file == "none") throw std::logic_error("The required model is not provided.\n");
    std::shared_ptr<float2DReg> mod = VECEXT::sepRead(mod_file);

// read in the 2D PEF parameters
    int pef_nx, pef_nz, pef_lead, winx, winz, version, pef_niter;
	float epsilon, weightx;
    MB::readParam(argc, argv, "n1", pef_nz, 6);
    MB::readParam(argc, argv, "n2", pef_nx, 3);
    MB::readParam(argc, argv, "lead", pef_lead,0);
	MB::readParam(argc, argv, "epsilon", epsilon,0.0);
	MB::readParam(argc, argv, "version", version,0);
	MB::readParam(argc, argv, "winx", winx,20);
	MB::readParam(argc, argv, "winz", winz,20);
	MB::readParam(argc, argv, "pef_niter", pef_niter,20);
    MB::readParam(argc, argv, "weightx", weightx,0);

// read in the initial PEF if provided
    std::string pef_file;
	std::shared_ptr<float2DReg> pef;
	MB::readParam(argc, argv, "pef", pef_file, "none");
	if (pef_file != "none"){
    	std::shared_ptr<float2DReg> pef0 = VECEXT::sepRead(pef_file);
		pef = pef0;
	}
	else{
		std::shared_ptr<float2DReg> pef0 (new float2DReg(pef_nz,pef_nx));
		pef = pef0;
	}

// read in other parameters
    int niter;
    float threshold;
    std::string objFunc, grad, res, solver;
    MB::readParam(argc, argv, "niter", niter, 10);
    MB::readParam(argc, argv, "threshold", threshold, 0.01);
    MB::readParam(argc, argv, "obj_func", objFunc,"none");
    MB::readParam(argc, argv, "grad", grad,"none");
    MB::readParam(argc, argv, "res", res,"none");
    MB::readParam(argc, argv, "solver", solver,"cgls");

// set the data vector
    std::shared_ptr<float2DReg> zero = mod->clone();
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

//set the operator
    std::shared_ptr<nsPef2D> Ns;
    switch (version){
		case 0:
        {
			std::shared_ptr<nsPef2D_v0> nsP( new nsPef2D_v0(pef, dat0, mod->getHyper(), dat->getHyper(), pef_lead, epsilon));
			Ns = nsP;
        }
			break;

		case 1:
        {
			std::shared_ptr<nsPef2D_v1> nsP( new nsPef2D_v1(pef, dat0, mod->getHyper(), dat->getHyper(), pef_lead, epsilon));
            Ns = nsP;
        }
			break;

		case 2:
        {
			std::shared_ptr<nsPef2D_v2> nsP( new nsPef2D_v2(pef, dat0, mod->getHyper(), dat->getHyper(), pef_lead, epsilon));
            nsP->setWeightx(weightx);
            Ns = nsP;
        }
			break;

        case 3:
        {
			std::shared_ptr<nsPef2D_v3> nsP( new nsPef2D_v3(pef, dat0, mod->getHyper(), dat->getHyper(), pef_lead, epsilon, weightx));
            Ns = nsP;
        }
			break;

        case 4:
        {
			std::shared_ptr<nsPef2D_v4> nsP( new nsPef2D_v4(pef, dat0, mod->getHyper(), dat->getHyper(), pef_lead, epsilon));
            Ns = nsP;
        }
			break;

        case 5:
        {
			std::shared_ptr<nsPef2D_v5> nsP( new nsPef2D_v5(pef, dat0, mod->getHyper(), dat->getHyper(), pef_lead, epsilon));
            nsP->setWeightx(weightx);
            Ns = nsP;
        }
			break;

        case 9:
        {
			std::shared_ptr<nsPef2D_v9> nsP( new nsPef2D_v9(pef, dat0, mod->getHyper(), dat->getHyper(), pef_lead, winx, winz, pef_niter));
			Ns = nsP;
        }
			break;

		default:
			throw std::logic_error("The selected version is not implemented");
    }

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
    sol->run(Ns.get(), mod, dat);
        
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