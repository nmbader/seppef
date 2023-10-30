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
#include "fkTransform.h"
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

// compute FK spectrum of the data
    axis X = dat->getHyper()->getAxis(2);
    axis Z = dat->getHyper()->getAxis(1);
    fkTransform fk(dat->getHyper());
    X.d = fk.dk();
    X.o = fk.ok();
    Z.d = fk.df();
    Z.o = 0;
    Z.n = Z.n /2 +1;
    std::shared_ptr<complex2DReg> vec (new complex2DReg(Z, X));
    fk.forward(false,dat,vec);
    std::shared_ptr<float2DReg> temp (new float2DReg(vec->getHyper()));
    std::shared_ptr<float2DReg> datfk (new float2DReg(vec->getHyper()));
    temp = CVECEXT::module(vec);

// revert the K axis
    for (int ix=0; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            (*datfk->_mat)[X.n-ix-1][iz] = (*temp->_mat)[ix][iz];
        }
    }


// write the data fk spectrum to disk
    VECEXT::sepWrite(datfk, output);
	
    return 0;
}