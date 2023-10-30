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
#include "seplib.h"

using namespace SEP;


// ========================================
// ================ Main ==================
// ========================================


int main(int argc, char **argv) {
    
    std::clog << "Starting program generate_random_noise \n";
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the parameters
    int nx, nz, xrand, zrand, seed;
    float ox, oz, dx, dz, min, max;
    std::string xlabel, zlabel;
    MB::readParam(argc, argv, "n1", nz, 1);
    MB::readParam(argc, argv, "n2", nx, 1);
    MB::readParam(argc, argv, "d1", dz, 1);
    MB::readParam(argc, argv, "d2", dx, 1);
    MB::readParam(argc, argv, "o1", oz, 0);
    MB::readParam(argc, argv, "o2", ox, 0);
    MB::readParam(argc, argv, "xlabel", xlabel, "Offset (m)");
    MB::readParam(argc, argv, "zlabel", zlabel, "Time (s)");
    MB::readParam(argc, argv, "xrand_inc", xrand, 0);
    MB::readParam(argc, argv, "zrand_inc", zrand, 0);
    MB::readParam(argc, argv, "min", min, -1);
    MB::readParam(argc, argv, "max", max, 1);
    MB::readParam(argc, argv, "seed", seed, 1987);

    srand(seed);

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
	
// filling data with random numbers, at random locations if xrand or zrand != 0
    for (int ix = 0; ix < nx; ix+=std::max(1, xrand*(rand()%10000)/10000)){
        for (int iz=0; iz < nz; iz+=std::max(1, zrand*(rand()%10000)/10000)){
            (*dat->_mat)[ix][iz] = min + (max-min)*(rand() % 10000)/10000.0;
        }
    }    

// write the data to disk
    VECEXT::sepWrite(dat, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program generate_random_noise; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}