#include <stdio.h>
#include <string>
#include <math.h>

#include "floatHyperExt.h"
#include "float2DReg.h"
#include "functions.h"
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program generate_monochromatic_wave \n";
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the parameters
    int nx, nz;
    float ox, oz, dx, dz, vel, amp, w, phase;
    std::string xlabel, zlabel;
    MB::readParam(argc, argv, "n1", nz, 1);
    MB::readParam(argc, argv, "n2", nx, 1);
    MB::readParam(argc, argv, "d1", dz, 1);
    MB::readParam(argc, argv, "d2", dx, 1);
    MB::readParam(argc, argv, "o1", oz, 0);
    MB::readParam(argc, argv, "o2", ox, 0);
    MB::readParam(argc, argv, "xlabel", xlabel, "X (m)");
    MB::readParam(argc, argv, "zlabel", zlabel, "Time (s)");
    MB::readParam(argc, argv, "vel", vel, 1500);
    MB::readParam(argc, argv, "amp", amp, 1);
    MB::readParam(argc, argv, "w", w, 0.1);
    MB::readParam(argc, argv, "phase", phase, 0.0);

    assert((w>=0) && (w<=1));


// initializing the data
	axis X, Z;
	X.o=ox; // x origin (in m)
	Z.o=oz; // z origin (in m)
	X.d=dx; // x sampling (in m)
	Z.d=dz; // (in sec)
	X.n=nx; // number of samples in x
	Z.n=nz; // number of samples in z
    X.label = xlabel;
    Z.label = zlabel;
   
	std::shared_ptr<float2DReg> dat (new float2DReg(Z, X));
	
// filling data with a monochromatic plane wave
    double arg;
    for (int ix=0; ix<X.n; ix++){
        for (int iz=0; iz<Z.n; iz++){
            arg = w*M_PI/Z.d*(iz*Z.d+Z.o - (ix*X.d+X.o)/vel);
            (*dat->_mat)[ix][iz] = amp*cos(arg+phase);
        }
    }
	

// write the data to disk
    VECEXT::sepWrite(dat, output);
	
    std::clog << "Ending program generate_monochromatic_wave\n";
    return 0;
}