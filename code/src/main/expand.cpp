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
    
    std::clog << "Starting program expand \n";
	srand ( time(0) ); 
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in the dimensions of the expansion
    bool xExpand, zExpand;
    MB::readParam(argc, argv, "n1", zExpand, true);
    MB::readParam(argc, argv, "n2", xExpand, true);

// expand the data
    std::shared_ptr<float2DReg> dat_exp;
    dat_exp = VECEXT::expand(dat,xExpand,zExpand);

// write the expanded data to disk
    VECEXT::sepWrite(dat_exp, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program expand; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}