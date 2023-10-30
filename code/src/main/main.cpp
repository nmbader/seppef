//
//  main.cpp
//  PEF
//
//  Created by Milad Bader on 11/13/18.
//  Copyright © 2018 Milad Bader. All rights reserved.
//

#include <time.h>
#include <chrono>
//#include <boost/tuple/tuple.hpp>
//#include <boost/foreach.hpp>

#include "gnuplot_iostream.h"

#include <sys/stat.h>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include "float1DReg.h"
#include "float2DReg.h"
#include "floatHyperExt.h"
#include "functions.h"
#include "convolution2D.h"
#include "chainOper2D.h"
#include "fxTransform.h"
#include "fkTransform.h"
#include "cgls.h"
#include "cglsReg.h"
#include "plotting.h"
#include "axis.h"
#include "hypercube.h"
#include "cglsSuper.h"
#include "pef2D.h"
#include "scaleAddOper2D.h"
#include "nsPef2D.h"
#include "expNsPef2D.h"
#include "sdls.h"
#include "varDatConv2D.h"
#include "varConv2D.h"
#include "expVarConv2D.h"
#include "varPef2D.h"
//#include "sigpack.h"
#include "doubleHyperExt.h"

#include "seplib.h"
//#include "sep3d.h"
//#include "sep_pars_external.h"
//#include "sepaux.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {

	initpar(argc,argv);

	std::string input, outmod, outdat;
	MB::readParam(argc, argv, "input", input,"in");
	MB::readParam(argc, argv, "outmod", outmod,"none");
	MB::readParam(argc, argv, "outdat", outdat,"none");

	srand(1234);

    std::shared_ptr<float2DReg> fil (new float2DReg(10,10));
	VECEXT::random(fil,-0.1,0.1);

	std::shared_ptr<float2DReg> dat (new float2DReg(1000,1000));
	VECEXT::random(dat,-0.1,0.1);
	varDatConv2D conv(dat, fil->getHyper(), 20, 20, 5);
	std::shared_ptr<float2DReg> mod (new float2DReg(conv.getDomain()));
	VECEXT::random(mod,-0.1,0.1);
	std::shared_ptr<float2DReg> temp = mod->clone();

	auto start = std::chrono::high_resolution_clock::now();
	for (int i=0; i<20; i++){
		mod=temp->clone();
		conv.forward(false, mod, dat);
		conv.adjoint(false, mod, dat);
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program; execution time (sec): " <<  duration.count()/1000000.0 << "\n";

	conv.dotProduct();

std::clog << dat->norm(2) << "\t" << mod->norm(2) << "\n";

	if (outdat != "none")
		VECEXT::sepWrite(dat, outdat);
	if (outmod != "none")	
		VECEXT::sepWrite(mod, outmod);

	/* cgls cg(1, 0.001);
    varPef2D p(dat->getHyper(), 10, 10, 20, 20, 5);
    p.estimate(dat, cg); */

	//VECEXT::iirBpButterworth(dat,0.05, 0.2,8);

	/* {
    plot1D plt(dat->getHyper());
    plt.plot(dat);
    } */

	/* sp::FIR_filt<double, double, double> fir_filt;
	arma::vec b;
	b = sp::fir1(100, 0.1);
  	fir_filt.set_coeffs(b);

	for (int i = 0; i< Z.n; i++){
		dat->getVals()[i] = fir_filt(dat->getVals()[i]);
	}  

	{
    plot1D plt(dat->getHyper());
    plt.plot(dat);
    } */



    return 0;
}