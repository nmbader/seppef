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

// Convert f-k spectrum to frequency-phase velocity spectrum using linear interpolation

int main(int argc, char **argv) {
    
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in parameters
    float vmax, dv;
    MB::readParam(argc, argv, "vmax", vmax, 1000);
    MB::readParam(argc, argv, "dv", dv, 10);

// compute FV spectrum of the data
    axis X1 = dat->getHyper()->getAxis(2);
    axis Z1 = dat->getHyper()->getAxis(1);
    axis X2 = Z1;
    axis Z2 = X1;
    
    Z2.n = 2*round(vmax/dv) + 1;
    vmax = round(vmax/dv)*dv;
    Z2.d = dv;
    Z2.o = -vmax;
    
    std::shared_ptr<float2DReg> vec (new float2DReg(Z2, X2));
    vec->zero();
    float kmax = abs(X1.o);
    float k, v, f, vmin, dk1, dk2;
    int ivmin, ik1, ik2;
    for (int ix=1; ix<X2.n; ix++){
        f = ix*X2.d;
        vmin = 2*M_PI*f/kmax;
        ivmin = ceil(vmin/dv);
        for (int iz=(Z2.n-1)/2+ivmin; iz<Z2.n; iz++){
            
            v = iz*Z2.d + Z2.o;
            k = 2*M_PI*f/v;
            ik1 = ceil((k-X1.o)/X1.d);
            ik2 = ik1 - 1;
            dk1 = ik1*X1.d+X1.o-k;
            dk2 = k-(ik2*X1.d+X1.o);
            (*vec->_mat)[ix][iz] = dk1/X1.d*(*dat->_mat)[ik1][ix]+dk2/X1.d*(*dat->_mat)[ik2][ix];
            (*vec->_mat)[ix][Z2.n-1-iz] = dk1/X1.d*(*dat->_mat)[X1.n-1-ik1][ix]+dk2/X1.d*(*dat->_mat)[X1.n-1-ik2][ix];
        }
    }


// write the data fk spectrum to disk
    VECEXT::sepWrite(vec, output);
	
    return 0;
}