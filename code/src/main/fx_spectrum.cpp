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
#include "fxTransform.h"
#include "complexHyperExt.h"
#include "complex2DReg.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();


// compute FX spectrum of the data
    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);
    fxTransform fx(dat->getHyper());
    Z.d = fx.df();
    Z.o = 0;
    Z.n = Z.n /2 +1;
    std::shared_ptr<complex2DReg> vec (new complex2DReg(Z, X));
    fx.forward(false,dat,vec);
    std::shared_ptr<float2DReg> datfx (new float2DReg(vec->getHyper()));
    datfx = CVECEXT::module(vec);
   
   
    VECEXT::sepWrite(datfx, output);
	
    return 0;
}