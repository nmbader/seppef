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
    
    std::clog << "Starting program cut \n";
	auto start = std::chrono::high_resolution_clock::now();
    initpar(argc,argv);
	std::string output;
	MB::readParam(argc, argv, "output", output,"out");

// read in the input data
    std::shared_ptr<float2DReg> dat = VECEXT::sepRead();

// read in the cut parameters
    int cut_n1_start, cut_n2_start, cut_n1_end, cut_n2_end;
    MB::readParam(argc, argv, "cut_n1_start", cut_n1_start, 0);
    MB::readParam(argc, argv, "cut_n2_start", cut_n2_start, 0);
    MB::readParam(argc, argv, "cut_n1_end", cut_n1_end, 0);
    MB::readParam(argc, argv, "cut_n2_end", cut_n2_end, 0);

// cut the data
    std::shared_ptr<float2DReg> mod = VECEXT::cut(dat, cut_n2_start, cut_n2_end, cut_n1_start, cut_n1_end);

// write the cut data to disk
    VECEXT::sepWrite(mod, output);
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::clog << "Ending program cut; execution time (sec): " <<  duration.count()/1000000.0 << "\n";
    return 0;
}