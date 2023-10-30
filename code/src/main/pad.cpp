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
    
    std::clog << "Starting program pad \n";
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in the padding parameters
    int pad_n1_start, pad_n2_start, pad_n1_end, pad_n2_end;
    MB::readParam(argc, argv, "pad_n1_start", pad_n1_start, 0);
    MB::readParam(argc, argv, "pad_n2_start", pad_n2_start, 0);
    MB::readParam(argc, argv, "pad_n1_end", pad_n1_end, 0);
    MB::readParam(argc, argv, "pad_n2_end", pad_n2_end, 0);

// pad the data
    std::shared_ptr<float2DReg> mod = VECEXT::pad(dat, pad_n2_start, pad_n2_end, pad_n1_start, pad_n1_end);

// write the padded data to disk
    VECEXT::sepWrite(mod, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program pad; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}